"""
utility functions for gene mapper class
"""

import pathlib
import sqlite3
import tempfile

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.ncbi_ingestion as ncbi_ingestion
import mmc_gene_mapper.create_db.species_ingestion as species_ingestion
import mmc_gene_mapper.create_db.bkbit_ingestion as bkbit_ingestion
import mmc_gene_mapper.create_db.ortholog_ingestion as ortholog_ingestion


def create_mapper_database(
        db_path,
        download_manager,
        tmp_dir,
        force_download,
        data_file_spec=None):

    with sqlite3.connect(db_path) as conn:
        data_utils.create_data_tables(conn)
        metadata_utils.create_metadata_tables(conn)

    ncbi_ingestion.ingest_ncbi_data(
        db_path=db_path,
        download_manager=download_manager,
        clobber=False
    )

    species_ingestion.ingest_species_data(
        db_path=db_path,
        download_manager=download_manager,
        tmp_dir=tmp_dir,
        force_download=force_download
    )

    if data_file_spec is not None:
        for file_spec in data_file_spec:
            if file_spec['type'] == 'bkbit':
                bkbit_ingestion.ingest_bkbit_genes(
                    db_path=db_path,
                    bkbit_path=file_spec['path'],
                    citation_name=file_spec['name']
                )
            elif file_spec['type'] == 'hmba_orthologs':
                ortholog_ingestion.ingest_hmba_orthologs(
                    db_path=db_path,
                    hmba_file_path=file_spec['path'],
                    citation_name=file_spec['name'],
                    baseline_species=file_spec.get('baseline_species', 'human'),
                    clobber=False
                )
            else:
                raise RuntimeError(
                    f"cannot parse file of type {file_spec['type']}"
                )

    # only after all data has been ingested
    with sqlite3.connect(db_path) as conn:
        data_utils.create_data_indexes(conn)
