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


@pytest.mark.parametrize("strict", [True, False])
def test_get_species_taxon(
        species_db_fixture,
        strict):

    assert query_utils.get_species_taxon(
        db_path=species_db_fixture,
        species_name='house mouse',
        strict=strict) == 10090

    assert query_utils.get_species_taxon(
        db_path=species_db_fixture,
        species_name='Mus musculus',
        strict=strict) == 10090

    assert query_utils.get_species_taxon(
        db_path=species_db_fixture,
        species_name='human',
        strict=strict) == 9606

    with pytest.raises(ValueError, match='species match name'):
        query_utils.get_species_taxon(
            db_path=species_db_fixture,
            species_name='dragon',
            strict=strict
        )

    if strict:
        with pytest.raises(ValueError, match='no species match'):
            query_utils.get_species_taxon(
                db_path=species_db_fixture,
                species_name='goblin',
                strict=strict
            )
    else:
        assert query_utils.get_species_taxon(
            db_path=species_db_fixture,
            species_name='goblin',
            strict=strict
        ) is None


@pytest.mark.parametrize("strict", [True, False])
def test_get_species_taxon_from_cursor(
        species_db_fixture,
        strict):

    with sqlite3.connect(species_db_fixture) as conn:
        cursor = conn.cursor()

        assert query_utils._get_species_taxon(
            cursor=cursor,
            species_name='house mouse',
            strict=strict) == 10090

        assert query_utils._get_species_taxon(
            cursor=cursor,
            species_name='Mus musculus',
            strict=strict) == 10090

        assert query_utils._get_species_taxon(
            cursor=cursor,
            species_name='human',
            strict=strict) == 9606

        with pytest.raises(ValueError, match='species match name'):
            query_utils._get_species_taxon(
                cursor=cursor,
                species_name='dragon',
                strict=strict
            )

        if strict:
            with pytest.raises(ValueError, match='no species match'):
                query_utils._get_species_taxon(
                    cursor=cursor,
                    species_name='goblin',
                    strict=strict
                )
        else:
            assert query_utils._get_species_taxon(
                cursor=cursor,
                species_name='goblin',
                strict=strict
            ) is None
