import pytest

import json
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.create_db.bkbit_ingestion as bkbit_ingestion
import mmc_gene_mapper.create_db.species_ingestion as species_ingestion


@pytest.fixture(scope='function')
def pre_bkbit_database_fixture(
        tmp_dir_fixture,
        species_file_fixture,
        dummy_download_mgr_fixture):
    """
    Return the path to a database that is prepped for
    bkbit ingestion (basically: species and authorities
    have been ingested)
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='bkbit_ingest_db_',
        suffix='.db',
        delete=True
    )

    species_ingestion.ingest_species_data(
        db_path=db_path,
        download_manager=dummy_download_mgr_fixture(),
        force_download=True,
        tmp_dir=tmp_dir_fixture
    )

    with sqlite3.connect(db_path) as conn:
        metadata_utils.create_authority_table(conn.cursor())
        metadata_utils.insert_authority(
            conn=conn,
            name="JABB",
            strict=True
        )

    return db_path


def test_read_bkbit_data(
        pre_bkbit_database_fixture,
        bkbit_data_fixture0):
    """
    Test that function to read a bkbit jsonld file correctly
    parses and returns the data
    """

    bkbit_data = bkbit_ingestion.read_bkbit_data(
        bkbit_path=bkbit_data_fixture0,
        db_path=pre_bkbit_database_fixture
    )

    with open(bkbit_data_fixture0, 'rb') as src:
        raw_data = json.load(src)

    assert bkbit_data[2] == "J001-2025"
    assert bkbit_data[0] == {
        "biolink:OrganismTaxon": raw_data['@graph'][0],
        "bican:GenomeAssembly": raw_data['@graph'][1],
        "bican:GenomeAnnotation": raw_data['@graph'][2]
    }

    expected_tuples = []
    for ii, el in enumerate(raw_data['@graph'][3:]):
        this = (
            0,
            ii,
            999,
            el['symbol'],
            el['source_id']
        )
        expected_tuples.append(this)
        if el['symbol'] != el['name']:
            this = (
                0,
                ii,
                999,
                el['name'],
                el['source_id']
            )
            expected_tuples.append(this)

    assert expected_tuples == bkbit_data[1]


def test_ingest_bkbit_data(
        pre_bkbit_database_fixture,
        bkbit_data_fixture0):
    """
    Test function to ingest bkbit data into gene table
    """
    with sqlite3.connect(pre_bkbit_database_fixture) as conn:
        cursor = conn.cursor()
        metadata_utils.create_citation_table(cursor)
        data_utils.create_gene_table(cursor)

    bkbit_ingestion.ingest_bkbit_genes(
        db_path=pre_bkbit_database_fixture,
        bkbit_path=bkbit_data_fixture0
    )

    with open(bkbit_data_fixture0, 'rb') as src:
        raw_data = json.load(src)

    expected_tuples = []
    for ii, el in enumerate(raw_data['@graph'][3:]):
        this = (
            0,
            ii,
            999,
            el['symbol'],
            el['source_id'],
            0
        )
        expected_tuples.append(this)
        if el['symbol'] != el['name']:
            this = (
                0,
                ii,
                999,
                el['name'],
                el['source_id'],
                0
            )
            expected_tuples.append(this)

    with sqlite3.connect(pre_bkbit_database_fixture) as conn:
        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                authority,
                id,
                species_taxon,
                symbol,
                identifier,
                citation
            FROM
                gene
            """
        ).fetchall()

    assert len(expected_tuples) == len(actual)
    assert len(set(actual)) == len(actual)
    assert set(actual) == set(expected_tuples)


def test_reingest_bkbit_data(
        pre_bkbit_database_fixture,
        bkbit_data_fixture0):
    """
    Test that an error occurs if you try to ingest two bkbit
    files claiming to be the same citation
    """
    with sqlite3.connect(pre_bkbit_database_fixture) as conn:
        cursor = conn.cursor()
        metadata_utils.create_citation_table(cursor)
        data_utils.create_gene_table(cursor)

    bkbit_ingestion.ingest_bkbit_genes(
        db_path=pre_bkbit_database_fixture,
        bkbit_path=bkbit_data_fixture0
    )

    with pytest.raises(ValueError, match="citation J001-2025 already exists"):
        bkbit_ingestion.ingest_bkbit_genes(
            db_path=pre_bkbit_database_fixture,
            bkbit_path=bkbit_data_fixture0
        )
