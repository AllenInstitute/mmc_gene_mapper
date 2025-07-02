import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.query_db.query as query_utils


@pytest.fixture(scope='module')
def species_db_fixture(tmp_dir_fixture):
    """
    Path to database with a species table in it
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='species_query_test_',
        suffix='.db'
    )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE NCBI_species (
                name STRING,
                id INTEGER
            )
            """
        )

        values = [
            ('house mouse', 10090),
            ('Mus musculus', 10090),
            ('human', 9606),
            ('dragon', 9607),
            ('dragon', 999)
        ]

        cursor.executemany(
            """
            INSERT INTO NCBI_species (
                name,
                id
            )
            VALUES (?, ?)
            """,
            values
        )

    return db_path


def test_get_species_taxon(
        species_db_fixture):

    with sqlite3.connect(species_db_fixture) as conn:
        cursor = conn.cursor()

        assert query_utils._get_species_taxon(
            cursor=cursor,
            species_name='house mouse') == 10090

        assert query_utils._get_species_taxon(
            cursor=cursor,
            species_name='Mus musculus') == 10090

        assert query_utils._get_species_taxon(
            cursor=cursor,
            species_name='human') == 9606

        with pytest.raises(ValueError, match='species match name'):
            query_utils._get_species_taxon(
                cursor=cursor,
                species_name='dragon'
            )

        with pytest.raises(ValueError, match='no species match for'):
            query_utils._get_species_taxon(
                cursor=cursor,
                species_name='nothingburger'
            )


def test_get_species_name(
        species_db_fixture):

    with sqlite3.connect(species_db_fixture) as conn:
        cursor = conn.cursor()

        assert query_utils._get_species_name(
            cursor=cursor,
            species_taxon=10090) in ('house mouse', 'Mus Musculus')

        assert query_utils._get_species_name(
            cursor=cursor,
            species_taxon=9606) == 'human'

        with pytest.raises(ValueError, match='species match for'):
            query_utils._get_species_name(
                cursor=cursor,
                species_taxon=123545
            )
