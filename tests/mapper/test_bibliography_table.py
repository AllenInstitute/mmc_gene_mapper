"""
These tests of the bibliography table assume the full schema
(with the has_symbols) column
"""
import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils



@pytest.fixture
def gene_table_fixture(
        tmp_dir_fixture):
    """
    Return path to database file with gene table created
    and populated for test

    bibliography truth
    
    species, authority, citation, has_symbols
    0        0          0         1
    0        0          1         0
    0        1          1         1
    1        1          0         1
    1        1          1         0
    """

    values = [
        (0, 0, 0, 0, 'a'),
        (0, 0, 0, 1, None),
        (0, 0, 0, 2, 'b'),
        (0, 0, 1, 0, None),
        (0, 0, 1, 1, None),
        (0, 0, 1, 1, None),
        (0, 1, 1, 0, None),
        (0, 1, 1, 0, 'c'),
        (0, 1, 1, 0, None),
        (1, 1, 0, 0, 'x'),
        (1, 1, 0, 1, None),
        (1, 1, 0, 1, 'y'),
        (1, 1, 1, 0, None),
        (1, 1, 1, 1, None)
    ]

    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='bibliograph_table_test_',
        suffix='.db'
    )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        data_utils.create_gene_table(cursor)
        cursor.executemany(
            """
            INSERT INTO gene (
                species_taxon,
                authority,
                citation,
                id,
                symbol
            )
            VALUES  (?, ?, ?, ?, ?)
            """,
            values
        )


    return db_path


def test_has_symbols(gene_table_fixture):
    """
    Make sure that the has_symbols column is correctly populated
    """
    mapper_utils.create_bibliography_table(gene_table_fixture)
    
    with sqlite3.connect(gene_table_fixture) as conn:
        cursor = conn.cursor()
        actual = cursor.execute(
             """
             SELECT
                 species_taxon,
                 authority,
                 citation,
                 has_symbols
             FROM species_bibliography
             """
        ).fetchall()

    expected = set(
        [(0, 0, 0, 1),
         (0, 0, 1, 0),
         (0, 1, 1, 1),
         (1, 1, 0, 1),
         (1, 1, 1, 0)]
    )
    assert len(actual) == len(expected)
    assert set(actual) == expected
