import pytest

import json
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils
import mmc_gene_mapper.query_db.query as query_utils


@pytest.fixture(scope='session')
def gene_table_data_fixture():
    """
    A list of dicts representing the records
    to be inserted into the test gene table
    """
    return [
        {'authority': 0, 'citation': 1, 'species_taxon': 88,
         'id': 0, 'identifier': 'a', 'symbol': 'b'},
        {'authority': 0, 'citation': 1, 'species_taxon': 88,
         'id': 1, 'identifier': 'c', 'symbol': 'd'},
        {'authority': 1, 'citation': 2, 'species_taxon': 88,
         'id': 0, 'identifier': 'a', 'symbol': 'b'},
        {'authority': 0, 'citation': 2, 'species_taxon': 99,
         'id': 5, 'identifier': 'e', 'symbol': 'f'},
        {'authority': 1, 'citation': 1, 'species_taxon': 77,
         'id': 6, 'identifier': 'g', 'symbol': 'h'},
        {'authority': 1, 'citation': 2, 'species_taxon': 88,
         'id': 10, 'identifier': 'y', 'symbol': 'z'},
        {'authority': 0, 'citation': 3, 'species_taxon': 99,
         'id': 11, 'identifier': 'u', 'symbol': 'w'}
    ]


@pytest.fixture(scope='function')
def gene_table_fixture(
        tmp_dir_fixture,
        gene_table_data_fixture):
    """
    Return the path to a database table with a gene
    table and a bogus bibliography table so that we
    can test the functions to create the bibliography
    table.
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='gene_for_bibio_',
        suffix='.db'
    )

    citation_values = sorted(
        set(
            [row['citation'] for row in gene_table_data_fixture]
        )
    )

    authority_values = sorted(
        set(
            [row['authority'] for row in gene_table_data_fixture]
        )
    )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        data_utils.create_gene_table(cursor)

        values = [
            (row['authority'],
             row['citation'],
             row['species_taxon'],
             row['id'],
             row['identifier'],
             row['symbol'])
            for row in gene_table_data_fixture
        ]

        cursor.executemany(
            """
            INSERT INTO gene (
                authority,
                citation,
                species_taxon,
                id,
                identifier,
                symbol
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            values
        )

        # cursor create bogus bibliograph table
        cursor.execute(
            """
            CREATE TABLE species_bibliography (
                id INTEGER,
                name STRING
            )
            """
        )
        cursor.executemany(
            """
            INSERT INTO species_bibliography (
                id,
                name
            ) VALUES (?,?)
            """,
            [(0, 'a'), (1, 'b')]
        )

        cursor.execute(
            """
            CREATE TABLE authority (
                name STRING,
                id INTEGER
            )
            """
        )
        cursor.executemany(
            """
            INSERT INTO authority(
                name,
                id
            )
            VALUES (?, ?)
            """,
            [(f"A{ii}", ii) for ii in authority_values]
        )

        cursor.execute(
            """
            CREATE TABLE citation (
                name STRING,
                id INTEGER,
                metadata STRING
            )
            """
        )
        cursor.executemany(
            """
            INSERT INTO citation (
                name,
                id,
                metadata
            ) VALUES (?, ?, ?)
            """,
            [
                (f"C{ii}", ii, json.dumps({'meta': f'M{ii}'}))
                for ii in citation_values
            ]
        )

    return db_path


def test_create_bibliography(gene_table_fixture):

    # make sure test database starts out with the bogus
    # bibliography table
    with sqlite3.connect(gene_table_fixture) as conn:
        cursor = conn.cursor()
        raw = cursor.execute(
            """
            SELECT * FROM species_bibliography
            """
        ).fetchall()
        assert len(raw[0]) == 2

    mapper_utils.create_bibliography_table(
        db_path=gene_table_fixture
    )

    with sqlite3.connect(gene_table_fixture) as conn:
        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                authority,
                species_taxon,
                citation
            FROM species_bibliography
            """
        ).fetchall()
        expected = set([
            (0, 88, 1),
            (1, 88, 2),
            (0, 99, 2),
            (1, 77, 1),
            (0, 99, 3)
        ])

        assert len(actual) == len(set(actual))
        assert set(actual) == set(expected)


def test_get_citation_from_bibliography(gene_table_fixture):
    mapper_utils.create_bibliography_table(
        db_path=gene_table_fixture
    )

    with sqlite3.connect(gene_table_fixture) as conn:
        cursor = conn.cursor()
        actual = query_utils.get_citation_from_bibliography(
            cursor=cursor,
            authority_idx=0,
            species_taxon=88
        )
        assert actual == {
            "name": "C1",
            "idx": 1,
            "metadata": {"meta": "M1"}
        }

        actual = query_utils.get_citation_from_bibliography(
            cursor=cursor,
            authority_idx=1,
            species_taxon=77
        )
        assert actual == {
            "name": "C1",
            "idx": 1,
            "metadata": {"meta": "M1"}
        }

        actual = query_utils.get_citation_from_bibliography(
            cursor=cursor,
            authority_idx=1,
            species_taxon=88
        )
        assert actual == {
            "name": "C2",
            "idx": 2,
            "metadata": {"meta": "M2"}
        }

        with pytest.raises(ValueError, match="2 citations associated with"):
            query_utils.get_citation_from_bibliography(
                cursor=cursor,
                authority_idx=0,
                species_taxon=99
            )

        with pytest.raises(ValueError, match="0 citations associated with"):
            query_utils.get_citation_from_bibliography(
                cursor=cursor,
                authority_idx=1,
                species_taxon=99
            )

        with pytest.raises(ValueError, match="0 citations associated with"):
            query_utils.get_citation_from_bibliography(
                cursor=cursor,
                authority_idx=15,
                species_taxon=99
            )


def test_get_authority_and_citation(gene_table_fixture):
    mapper_utils.create_bibliography_table(
        db_path=gene_table_fixture
    )

    with sqlite3.connect(gene_table_fixture) as conn:
        actual = query_utils.get_authority_and_citation(
            conn=conn,
            authority_name="A0",
            species_taxon=88
        )
        assert actual == {
            "authority": {
                "name": "A0",
                "idx": 0
            },
            "citation": {
                "name": "C1",
                "idx": 1,
                "metadata": {"meta": "M1"}
            }
        }

        actual = query_utils.get_authority_and_citation(
            conn=conn,
            authority_name="A1",
            species_taxon=77
        )
        assert actual == {
            "authority": {
                "name": "A1",
                "idx": 1
            },
            "citation": {
                "name": "C1",
                "idx": 1,
                "metadata": {"meta": "M1"}
            }
        }

        actual = query_utils.get_authority_and_citation(
            conn=conn,
            authority_name="A1",
            species_taxon=88
        )
        assert actual == {
            "authority": {
                "name": "A1",
                "idx": 1
            },
            "citation": {
                "name": "C2",
                "idx": 2,
                "metadata": {"meta": "M2"}
            }
        }

        with pytest.raises(ValueError, match="2 citations associated with"):
            query_utils.get_authority_and_citation(
                conn=conn,
                authority_name="A0",
                species_taxon=99
            )

        with pytest.raises(ValueError, match="0 citations associated with"):
            query_utils.get_authority_and_citation(
                conn=conn,
                authority_name="A1",
                species_taxon=99
            )

        msg = "does not exist in this database"
        with pytest.raises(ValueError, match=msg):
            query_utils.get_authority_and_citation(
                conn=conn,
                authority_name="A15",
                species_taxon=99
            )
