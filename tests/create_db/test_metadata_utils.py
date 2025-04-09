"""
tests for create_db/metadata_tables.py
"""
import pytest

import json
import shutil
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils


@pytest.fixture(scope='module')
def citation_dataset_fixture():
    """
    Data to be ingested into test citation table
    """
    return [
        {"name": "A",
         "metadata": {"A": 1, "B": {"x": "a", "y": 2}},
         "idx": 0},
        {"name": "B",
         "metadata": {"some": "metadata"},
         "idx": 1},
        {"name": "C",
         "metadata": {"other": "metadata"},
         "idx": 2},
        {"name": "C",
         "metadata": {"uh": "oh"},
         "idx": 3}
    ]


@pytest.fixture(scope='module')
def citation_table_fixture(
        citation_dataset_fixture,
        tmp_dir_fixture):
    """
    Return path to sqlite database with citation
    table to test against.
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='citation_table_',
        suffix='.db'
    )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE citation (
                id INTEGER,
                name STRING,
                metadata STRING
            )
            """
        )
        values = [
            (row['idx'], row['name'], json.dumps(row['metadata']))
            for row in citation_dataset_fixture
        ]
        cursor.executemany(
            """
            INSERT INTO citation (
                id,
                name,
                metadata
            ) VALUES (?, ?, ?)
            """,
            values
        )

    return db_path


def test_get_citation(
        citation_dataset_fixture,
        citation_table_fixture):

    with sqlite3.connect(citation_table_fixture) as conn:
        a_result = metadata_utils.get_citation(
            conn=conn,
            name="A")
        assert a_result == citation_dataset_fixture[0]

        a_result = metadata_utils.get_citation(
            conn=conn,
            name="B")
        assert a_result == citation_dataset_fixture[1]

        for strict in (True, False):
            with pytest.raises(ValueError, match="More than one citation"):
                metadata_utils.get_citation(
                    conn=conn,
                    name="C",
                    strict=strict
                )

        assert metadata_utils.get_citation(
            conn=conn,
            name="D",
            strict=False) is None

        with pytest.raises(ValueError, match="No citation"):
            metadata_utils.get_citation(
                conn=conn,
                name="D",
                strict=True
            )
