"""
test species_ingestion.py using simulated data
"""
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.species_ingestion as species_ingestion


def test_species_ingestion_smoke(
        tmp_dir_fixture,
        species_file_fixture,
        species_id_fixture):
    """
    Just try ingesting the species data and see what comes out
    """

    # define dummy downloan manager class
    class DummyDownloadManager(object):

        def get_file(
                self,
                host,
                src_path,
                force_download):
            if src_path.endswith('new_taxdump.tar.gz'):
                return {
                    'local_path': species_file_fixture['tar']
                }
            elif src_path.endswith('new_taxdump.tar.gz.md5'):
                return {
                    'local_path': species_file_fixture['hash']
                }
            else:
                raise RuntimeError(
                    "MockDownloadManager cannot handle src_path "
                    f"{src_path}"
                )

    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='species_db_',
        suffix='.db',
        delete=True
    )

    species_ingestion.ingest_species_data(
        db_path=db_path,
        download_manager=DummyDownloadManager(),
        force_download=True,
        tmp_dir=tmp_dir_fixture)

    expected_pairs = set([
        (key, species_id_fixture[key])
        for key in species_id_fixture
    ])

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        actual_pairs = set(
            cursor.execute(
                """
                SELECT name, id FROM NCBI_species
                """
            ).fetchall()
        )

    assert actual_pairs == expected_pairs
