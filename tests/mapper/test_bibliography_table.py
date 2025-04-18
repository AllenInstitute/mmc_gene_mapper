"""
These tests of the bibliography table assume the full schema
(with the has_symbols) column
"""
import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils
import mmc_gene_mapper.query_db.query as query_utils


@pytest.fixture
def gene_table_fixture(
        tmp_dir_fixture):
    """
    Return path to database file with gene table created
    and populated for test

    bibliography truth

    species, authority, citation, has_symbols
    0        0          0         0
    0        0          1         1
    0        1          1         1
    1        1          0         1
    1        1          1         0
    1        2          0         0
    2        3          2         1
    2        3          3         1
    """

    values = [
        (0, 0, 1, 0, 'a'),
        (0, 0, 1, 1, None),
        (0, 0, 1, 2, 'b'),
        (0, 0, 0, 0, None),
        (0, 0, 0, 1, None),
        (0, 0, 0, 1, None),
        (0, 1, 1, 0, None),
        (0, 1, 1, 0, 'c'),
        (0, 1, 1, 0, None),
        (1, 1, 0, 0, 'x'),
        (1, 1, 0, 1, None),
        (1, 1, 0, 1, 'y'),
        (1, 1, 1, 0, None),
        (1, 1, 1, 1, None),
        (1, 2, 0, 1, None),
        (1, 2, 0, 2, None),
        (2, 3, 2, 0, 'a'),
        (2, 3, 3,  0, 'b')
    ]

    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='bibliograph_table_test_',
        suffix='.db'
    )
    with sqlite3.connect(db_path) as conn:
        metadata_utils.create_metadata_tables(conn)
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

        for ii in range(4):
            metadata_utils.insert_citation(
                conn=conn,
                name=f'CITATION_{ii}',
                metadata_dict={'a': ['field', ii]}
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
        [(0, 0, 0, 0),
         (0, 0, 1, 1),
         (0, 1, 1, 1),
         (1, 1, 0, 1),
         (1, 1, 1, 0),
         (1, 2, 0, 0),
         (2, 3, 2, 1),
         (2, 3, 3, 1)]
    )
    assert len(actual) == len(expected)
    assert set(actual) == expected


def test_requires_symbols(gene_table_fixture):
    """
    Test that query_utils.get_citation_from_bibliography
    behaves correctly for different values of require_symbols
    """
    mapper_utils.create_bibliography_table(gene_table_fixture)
    error_msg = 'citations associated with authority'

    with sqlite3.connect(gene_table_fixture) as conn:
        cursor = conn.cursor()

        actual = query_utils.get_citation_from_bibliography(
            cursor=cursor,
            species_taxon=1,
            authority_idx=2,
            require_symbols=False
        )
        expected = {
            'idx': 0,
            'name': 'CITATION_0',
            'metadata': {'a': ['field', 0]}
        }
        assert actual == expected

        # there are no citations with symbosl for
        # species=1, authority=2
        with pytest.raises(ValueError, match=error_msg):
            query_utils.get_citation_from_bibliography(
                cursor=cursor,
                species_taxon=1,
                authority_idx=2,
                require_symbols=True
            )

        # check that if require_symbols=False, species=0, authority=0
        # returns the citation with symbols
        for require_symbols in (True, False):
            actual = query_utils.get_citation_from_bibliography(
                cursor=cursor,
                species_taxon=0,
                authority_idx=0,
                require_symbols=require_symbols
            )
            expected = {
                'idx': 1,
                'name': 'CITATION_1',
                'metadata': {'a': ['field', 1]}
            }
            assert actual == expected
