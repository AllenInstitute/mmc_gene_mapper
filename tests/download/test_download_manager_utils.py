import pytest

import pathlib
import sqlite3
import tempfile
import unittest.mock

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.download.download_manager_utils as mgr_utils


@pytest.fixture(scope='session')
def tmp_dir_fixture(
        tmp_path_factory):
    result = tmp_path_factory.mktemp(
        'mmc_gene_mapper_download_manager_'
    )
    yield result
    file_utils.clean_up(result)


def test_create_download_db(tmp_dir_fixture):
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_fixture)
    existing_file = file_utils.mkstemp_clean(
        dir=tmp_dir,
        prefix='existing_db_',
        suffix='.db',
        delete=False
    )

    with pytest.raises(ValueError, match="already exists"):
        mgr_utils.create_download_db(existing_file)

    valid_path = pathlib.Path(tmp_dir) / 'not_a_file.db'
    assert not valid_path.exists()
    mgr_utils.create_download_db(valid_path)
    assert valid_path.is_file()


@pytest.mark.parametrize(
    "file_name, host, src, expected",
    [("remove_record_db0.db",
      "h0",
      "s0",
      [("h1", "s0", "d/e/f", "ghijk", "lmnop"),
       ("h1", "s2", "j/k/l", "mnopq", "rstuv")
       ]
      ),
     ("remove_record_db1.db",
      "h1",
      "s0",
      [("h0", "s0", "a/b/c", "efghi", "jklmn"),
       ("h0", "s0", "g/h/i", "jklmn", "opqrs"),
       ("h1", "s2", "j/k/l", "mnopq", "rstuv")
       ]
      ),
     ]
)
def test_remove_records(
        file_name,
        host,
        src,
        expected,
        tmp_dir_fixture):
    """
    Test that mgr_utils.remove_reords really does remove all
    the records matching the input fields
    """

    db_path = pathlib.Path(tmp_dir_fixture) / file_name
    mgr_utils.create_download_db(db_path)

    values = [
        ("h0", "s0", "a/b/c", "efghi", "jklmn"),
        ("h1", "s0", "d/e/f", "ghijk", "lmnop"),
        ("h0", "s0", "g/h/i", "jklmn", "opqrs"),
        ("h1", "s2", "j/k/l", "mnopq", "rstuv")
    ]

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.executemany(
            """
            INSERT INTO downloads (
                host,
                src_path,
                local_path,
                hash,
                downloaded_on
            ) VALUES (?, ?, ?, ?, ?)
            """,
            values
        )

    mgr_utils.remove_record(
        db_path=db_path,
        host=host,
        src_path=src
    )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        result = cursor.execute(
            """
            SELECT
                host,
                src_path,
                local_path,
                hash,
                downloaded_on
            FROM
                downloads
            """
        ).fetchall()
    assert set(result) == set(expected)


def test_insert_record(tmp_dir_fixture):

    db_path = pathlib.Path(tmp_dir_fixture) / "insert_record.db"
    local_file = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='data_file_',
        suffix='.txt',
        delete=False
    )
    local_hash = file_utils.hash_from_path(local_file)
    mgr_utils.create_download_db(db_path)

    values = [
        ("h0", "s0", "a/b/c", "efghi", "jklmn"),
        ("h1", "s0", "d/e/f", "ghijk", "lmnop"),
        ("h0", "s0", "g/h/i", "jklmn", "opqrs"),
        ("h1", "s2", "j/k/l", "mnopq", "rstuv")
    ]

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.executemany(
            """
            INSERT INTO downloads (
                host,
                src_path,
                local_path,
                hash,
                downloaded_on
            ) VALUES (?, ?, ?, ?, ?)
            """,
            values
        )

    to_replace = "mmc_gene_mapper.utils.timestamp.get_timestamp"

    def dummy_timestamp():
        return "test-timestamp"

    with unittest.mock.patch(to_replace, dummy_timestamp):
        mgr_utils.insert_record(
            db_path=db_path,
            host="h2",
            src_path="s3",
            local_path=local_file
        )

    expected = values + [
        ("h2", "s3", local_file, local_hash, "test-timestamp")
    ]
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                host,
                src_path,
                local_path,
                hash,
                downloaded_on
            FROM
                downloads
            """
        ).fetchall()

    assert set(actual) == set(expected)


def test_insert_record_failures(tmp_dir_fixture):
    """
    Test that correct exception is raised if you try to
    insert a record in which local_path does not point to a valid
    file
    """
    db_path = pathlib.Path(tmp_dir_fixture) / "insert_record_bad.db"
    mgr_utils.create_download_db(db_path)

    msg = "mmc_gene_mapper/dummy/file.txt is not a file"
    with pytest.raises(file_utils.NotAFileError, match=msg):
        mgr_utils.insert_record(
            db_path=db_path,
            host='h5',
            src_path='s8',
            local_path="mmc_gene_mapper/dummy/file.txt"
        )


def test_mgr_utils_get_records(
        tmp_dir_fixture):

    db_path = pathlib.Path(tmp_dir_fixture) / "mgr_utils_get_record_test.db"
    mgr_utils.create_download_db(db_path)

    values = [
        ("h0", "s0", "a/b/c", "efghi", "jklmn"),
        ("h1", "s0", "d/e/f", "ghijk", "lmnop"),
        ("h0", "s0", "g/h/i", "jklmn", "opqrs"),
        ("h1", "s2", "j/k/l", "mnopq", "rstuv")
    ]

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.executemany(
            """
            INSERT INTO downloads (
                host,
                src_path,
                local_path,
                hash,
                downloaded_on
            ) VALUES (?, ?, ?, ?, ?)
            """,
            values
        )

    assert mgr_utils.get_record(
        db_path=db_path,
        host="h4",
        src_path="src4"
    ) == []

    actual = mgr_utils.get_record(
        db_path=db_path,
        host="h0",
        src_path="s0"
    )

    expected = [
        {"host": "h0",
         "src_path": "s0",
         "local_path": "a/b/c",
         "hash": "efghi",
         "downloaded_on": "jklmn"},
        {"host": "h0",
         "src_path": "s0",
         "local_path": "g/h/i",
         "hash": "jklmn",
         "downloaded_on": "opqrs"}
    ]

    assert len(actual) == len(expected)
    assert expected[0] in actual
    assert expected[1] in actual

    actual = mgr_utils.get_record(
        db_path=db_path,
        host="h1",
        src_path="s0"
    )

    expected = [
        {"host": "h1",
         "src_path": "s0",
         "local_path": "d/e/f",
         "hash": "ghijk",
         "downloaded_on": "lmnop"}
    ]

    assert actual == expected
