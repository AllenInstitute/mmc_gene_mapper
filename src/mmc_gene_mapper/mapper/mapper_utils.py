"""
utility functions for gene mapper class
"""

import numpy as np
import sqlite3

import mmc_gene_mapper.create_db.utils as db_utils
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
                    bkbit_path=file_spec['path']
                )

    if data_file_spec is not None:
        for file_spec in data_file_spec:
            if file_spec['type'] == 'bkbit':
                continue
            elif file_spec['type'] == 'hmba_orthologs':
                ortholog_ingestion.ingest_hmba_orthologs(
                    db_path=db_path,
                    hmba_file_path=file_spec['path'],
                    citation_name=file_spec['name'],
                    clobber=False
                )
            else:
                raise RuntimeError(
                    f"cannot parse file of type {file_spec['type']}"
                )

    # only after all data has been ingested
    with sqlite3.connect(db_path) as conn:
        data_utils.create_data_indexes(conn)


def create_bibliography_table(
        db_path):
    """
    Create species_bibliography which lists all
    combinations of authority, citation, species in the
    gene table.
    """
    print("=======CREATING BIBLIOGRAPHY TABLE=======")
    table_name = "species_bibliography"
    index_name = "species_bibliography_idx"
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        db_utils.delete_index(
            cursor=cursor,
            idx_name=index_name
        )
        cursor.execute(
            f"""
            DROP TABLE IF EXISTS {table_name}
            """
        )
        cursor.execute(
            f"""
            CREATE TABLE {table_name} (
                authority INTEGER,
                citation INTEGER,
                species_taxon INTEGER,
                has_symbols INTEGER
            )
            """
        )
        raw = cursor.execute(
            """
            SELECT
                authority,
                citation,
                species_taxon
            FROM gene
            GROUP BY
                authority,
                citation,
                species_taxon
            """
        ).fetchall()

        has_symbols = cursor.execute(
            """
            SELECT
                authority,
                citation,
                species_taxon,
                COUNT(symbol)
            FROM gene
            WHERE
                symbol IS NOT NULL
            GROUP BY
                authority,
                citation,
                species_taxon
            """
        ).fetchall()

        has_symbols = {
            row[:-1]: row[-1]
            for row in has_symbols
        }

        values = [
            (*row, 1)
            if row in has_symbols and has_symbols[row] > 0
            else (*row, 0)
            for row in raw
        ]

        cursor.executemany(
            f"""
            INSERT INTO {table_name} (
                authority,
                citation,
                species_taxon,
                has_symbols
            )
            VALUES (?, ?, ?, ?)
            """,
            values
        )
        db_utils.create_index(
            cursor=cursor,
            idx_name=index_name,
            table_name=table_name,
            column_tuple=("species_taxon", "authority", "has_symbols")
        )


def apply_mapping(
        gene_list,
        mapping):
    """
    Apply a mapping to a gene list, only keeping 1:1 matches.

    Parameters
    ----------
    gene_list:
        list of str. The input genes
    mapping:
        a dict mapping the genes in gene_list
        to the lists of matches produced by the MMCGeneMapper

    Returns
    -------
    A dict
        'failure_log': tally of how many genes failed and why
        'gene_list': the mapped gene list

    Notes
    -----
    Genes that fail to map in a 1:1 way will be assigned
    a name like UNMAPPABLE_{REASON}_{ii}

    Genes that are degenerate (i.e. that map to the same
    identifier) will also be marked as unmappable
    """
    failure_log = {
        'zero matches': 0,
        'many matches': 0,
        'degenerate matches': 0
    }

    new_gene_list = []
    for gene in gene_list:
        this = mapping[gene]
        if len(this) == 1:
            assn = this[0]
        else:
            if len(this) == 0:
                ct = failure_log['zero matches']
                assn = f'UNMAPPABLE_NO_MATCH_{ct}'
                failure_log['zero matches'] += 1
            else:
                ct = failure_log['many matches']
                assn = f'UNMAPPABLE_MANY_MATCHES_{ct}'
                failure_log['many matches'] += 1
        new_gene_list.append(assn)

    # make sure genes have unique names
    unq, ct = np.unique(new_gene_list, return_counts=True)
    degen = (ct > 1)
    unq = unq[degen]
    degen_gene_to_idx = {
        u: ii for ii, u in enumerate(unq)
    }
    degen_gene_to_ct = {
        g: 0 for g in degen_gene_to_idx
    }
    if len(degen_gene_to_idx) > 0:
        failure_log['degenerate matches'] = len(degen_gene_to_idx)*2
        for ii in range(len(new_gene_list)):
            gene = new_gene_list[ii]
            if gene in degen_gene_to_idx:
                pair = degen_gene_to_idx[gene]
                ct = degen_gene_to_ct[gene]
                assn = f'UNMAPPABLE_DEGENERATE_{pair}_{ct}'
                new_gene_list[ii] = assn
                degen_gene_to_ct[gene] += 1

    return {
        'failure_log': failure_log,
        'gene_list': new_gene_list
    }
