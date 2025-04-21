import pytest

import pathlib
import sqlite3
import tempfile

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
    [("db0.db",
      "h0",
      "s0",
      [("h1", "s0", "d/e/f", "ghijk", "lmnop"),
       ("h1", "s2", "j/k/l", "mnopq", "rstuv")
       ]
      ),
     ("db1.db",
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
