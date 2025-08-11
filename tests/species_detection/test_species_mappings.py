import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.mapper.species_detection as species_detection


@pytest.fixture(scope='session')
def species_mapping_db_fixture(tmp_dir_fixture):
    """
    Create a database with minimal tables to test the
    mapping of gene identifiers by species and authority
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix="species_mapping_",
        suffix=".db"
    )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE authority (
                id INTEGER,
                name TEXT
            )
            """
        )
        cursor.execute(
            """
            CREATE TABLE gene (
                identifier TEXT,
                symbol TEXT,
                species_taxon INTEGER,
                authority INTEGER,
                citation INTEGER
            )
            """
        )

        cursor.execute(
            """
            INSERT INTO authority (id, name)
            VALUES
                (0, 'NCBI'),
                (1, 'ENSEMBL')
            """
        )

        cursor.execute(
            """
            INSERT INTO gene (
                identifier,
                symbol,
                species_taxon,
                authority,
                citation)
            VALUES
                ('aaa', 'xxx', 0, 0, 1),
                ('aaa', 'zzz', 1, 0, 1),
                ('bbb', 'yyy', 0, 1, 1),
                ('ddd', 'xxx', 1, 1, 1),
                ('ddd', 'zzz', 0, 0, 1),
                ('aaa', 'xxx', 0, 0, 2)
            """
        )
    return db_path


def test_map_from_identifiers(species_mapping_db_fixture):

    actual = species_detection.species_from_identifier(
        db_path=species_mapping_db_fixture,
        gene_list=['aaa', 'bbb', 'ccc', 'xxx', 'eee', 'ddd'],
        chunk_size=2
    )

    assert len(actual) == 3
    assert len(actual['aaa']) == 2
    assert {'species_taxon': 0, 'authority': 'NCBI'} in actual['aaa']
    assert {'species_taxon': 1, 'authority': 'NCBI'} in actual['aaa']
    assert len(actual['bbb']) == 1
    assert {'species_taxon': 0, 'authority': 'ENSEMBL'} in actual['bbb']
    assert len(actual['ddd']) == 2
    assert {'species_taxon': 0, 'authority': 'NCBI'} in actual['ddd']
    assert {'species_taxon': 1, 'authority': 'ENSEMBL'} in actual['ddd']


def test_map_from_symbols(species_mapping_db_fixture):

    actual = species_detection.species_from_symbol(
        db_path=species_mapping_db_fixture,
        gene_list=['xxx', 'aaa', 'bbb', 'yyy',
                   'ccc', 'zzz', 'ddd', 'eee', 'fff'],
        chunk_size=2
    )

    assert len(actual) == 3
    assert len(actual['xxx']) == 2
    assert {'species_taxon': 0, 'authority': 'NCBI'} in actual['xxx']
    assert {'species_taxon': 1, 'authority': 'ENSEMBL'} in actual['xxx']
    assert len(actual['yyy']) == 1
    assert {'species_taxon': 0, 'authority': 'ENSEMBL'} in actual['yyy']
    assert len(actual['zzz']) == 2
    assert {'species_taxon': 1, 'authority': 'NCBI'} in actual['zzz']
    assert {'species_taxon': 0, 'authority': 'NCBI'} in actual['zzz']
