import gzip
import json
import pandas as pd
import pathlib
import sqlite3
import time

import mmc_gene_mapper
import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.utils.str_utils as str_utils
import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.create_db.utils as db_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.query_db.query as db_query
import mmc_gene_mapper.download.download_utils as download_utils


def ingest_ncbi_data(
        db_path,
        download_manager,
        clobber=True,
        force_download=False):

    metadata_dict = dict()
    host = 'ftp.ncbi.nlm.nih.gov'
    file_list = [
        'gene/DATA/gene_info.gz',
        'gene/DATA/gene2ensembl.gz',
        'gene/DATA/gene_orthologs.gz'
    ]

    for src_path in file_list:
        record = download_manager.get_file(
            host=host,
            src_path=src_path,
            force_download=force_download
        )
        fname = pathlib.Path(src_path).name
        metadata_dict[fname] = record


    ensembl_path = metadata_dict['gene2ensembl.gz'].pop('local_path')
    ortholog_path = metadata_dict['gene_orthologs.gz'].pop('local_path')
    gene_info_path = metadata_dict['gene_info.gz'].pop('local_path')

    _ingest_ncbi_data(
        db_path=db_path,
        gene_info_path=gene_info_path,
        ensembl_path=ensembl_path,
        ortholog_path=ortholog_path,
        metadata_dict=metadata_dict,
        clobber=clobber,
        citation_name='NCBI'
    )


def _ingest_ncbi_data(
        db_path,
        gene_info_path,
        ortholog_path,
        ensembl_path,
        metadata_dict,
        clobber=False,
        citation_name='NCBI'):

    file_utils.assert_is_file(ensembl_path)
    file_utils.assert_is_file(ortholog_path)
    file_utils.assert_is_file(gene_info_path)

    db_exists = False
    db_path = pathlib.Path(db_path)
    if db_utils.check_existence(db_path):
        db_exists = True

    with sqlite3.connect(db_path) as conn:
        if not db_exists:
            data_utils.create_data_tables(conn)
            metadata_utils.create_metadata_tables(conn)

        citation_idx = metadata_utils.insert_unique_citation(
            conn=conn,
            name=citation_name,
            metadata_dict=metadata_dict,
            clobber=clobber
        )

        auth_idx = metadata_utils.insert_unique_authority(
            conn=conn,
            name='NCBI'
        )

        metadata_utils.insert_authority(
            conn=conn,
            name='ENSEMBL'
        )

        ingest_gene_info(
            conn=conn,
            data_path=gene_info_path,
            authority_idx=auth_idx,
            citation_idx=citation_idx,
        )
        ingest_gene_to_ensembl(
            conn=conn,
            data_path=ensembl_path,
            citation_idx=citation_idx
        )
        ingest_orthologs(
            conn=conn,
            data_path=ortholog_path,
            citation_idx=citation_idx
        )


def ingest_gene_info(conn, data_path, authority_idx, citation_idx):
    t0 = time.time()
    print('=======INGESTING GENE INFO=======')
    data_path = pathlib.Path(data_path)
    cursor = conn.cursor()
    chunk_size = 5000000
    query = """
    INSERT INTO gene (
        authority,
        species_taxon,
        id,
        symbol,
        identifier,
        citation
    ) VALUES (?, ?, ?, ?, ?, ?)
    """
    i0 = 0

    #with open_fn(data_path, mode=mode) as src:
    with pd.read_csv(data_path,
                     delimiter='\t',
                     chunksize=chunk_size,
                     usecols=['#tax_id', 'GeneID', 'Symbol']) as src:
        for chunk in src:
            values = [
                (authority_idx,
                 int(tax_id),
                 int(gene_id),
                 gene_symbol,
                 f'NCBIGene:{gene_id}',
                 citation_idx)
                for tax_id, gene_id, gene_symbol in zip(
                        chunk['#tax_id'].values,
                        chunk['GeneID'].values,
                        chunk['Symbol'].values)
            ]

            if len(values) == 0:
                break

            print(f'    chunk {i0:.2e}, {i0+len(values):.2e}')

            i0 += len(chunk)
            cursor.executemany(query, values)
            conn.commit()
    dur = (time.time()-t0)/60.0
    print(f'=======INGESTING gene_info TOOK {dur:.2e} minutes=======')


def ingest_gene_to_ensembl(conn, data_path, citation_idx):
    t0 = time.time()
    print('=======INGESTING GENE TO ENSEMBL=======')

    ncbi_idx = metadata_utils.get_authority(
        conn=conn,
        name='NCBI'
    )["idx"]

    ensembl_idx = metadata_utils.get_authority(
        conn=conn,
        name='ENSEMBL'
    )["idx"]

    cursor = conn.cursor()
    data = pd.read_csv(data_path, delimiter='\t')
    n_rows = len(data)
    chunk_size = 100000

    query = """
    INSERT INTO gene_equivalence (
        species_taxon,
        authority0,
        gene0,
        authority1,
        gene1,
        citation
    )
    VALUES (?, ?, ?, ?, ?, ?)
    """
    uploaded_pairs = set()
    with pd.read_csv(
            data_path,
            delimiter='\t',
            chunksize=chunk_size,
            usecols=['#tax_id', 'GeneID', 'Ensembl_gene_identifier'],
            dtype={'#tax_id': int,
                   'GeneID': int,
                   'Ensembl_gene_identifier': str}) as src:

        for chunk in src:
            values = []
            for taxon_id, ncbi_gene_id, raw_ens_id in zip(
                    chunk['#tax_id'].values,
                    chunk['GeneID'].values,
                    chunk['Ensembl_gene_identifier'].values):

                ncbi_gene_id = int(ncbi_gene_id)
                taxon_id = int(taxon_id)

                try:
                    ensembl_gene_id = str_utils.int_from_identifier(
                        raw_ens_id
                    )
                except ValueError:
                    continue

                pair = (ncbi_gene_id, ensembl_gene_id)
                if pair in uploaded_pairs:
                    continue
                uploaded_pairs.add(pair)

                values.append(
                    (taxon_id,
                     ncbi_idx,
                     ncbi_gene_id,
                     ensembl_idx,
                     ensembl_gene_id,
                     citation_idx)
                )
                values.append(
                    (taxon_id,
                     ensembl_idx,
                     ensembl_gene_id,
                     ncbi_idx,
                     ncbi_gene_id,
                     citation_idx)
                )

            cursor.executemany(query, values)
            conn.commit()
    dur = (time.time()-t0)/60.0
    print(f'=======INGESTING gene2ensembl TOOK {dur:.2e} minutes=======')


def ingest_orthologs(
        conn,
        data_path,
        citation_idx):
    t0 = time.time()
    print('=======INGESTING ORTHOLOGS=======')

    ncbi_idx = metadata_utils.get_authority(
        conn=conn,
        name='NCBI'
    )["idx"]

    cursor = conn.cursor()
    data = pd.read_csv(data_path, delimiter='\t')
    data = data[data['relationship'] == 'Ortholog']
    n_rows = len(data)
    chunk_size = 50000

    query = """
    INSERT INTO gene_ortholog(
        authority,
        species0,
        gene0,
        species1,
        gene1,
        citation
    )
    VALUES (?, ?, ?, ?, ?, ?)
    """

    for i0 in range(0, n_rows, chunk_size):
        chunk = data.iloc[i0:i0+chunk_size].to_dict(orient='records')
        values = [
            (ncbi_idx,
             int(row['#tax_id']),
             int(row['GeneID']),
             int(row['Other_tax_id']),
             int(row['Other_GeneID']),
             citation_idx)
            for row in chunk
        ]
        values += [
            (ncbi_idx,
             int(row['Other_tax_id']),
             int(row['Other_GeneID']),
             int(row['#tax_id']),
             int(row['GeneID']),
             citation_idx)
            for row in chunk
        ]
        cursor.executemany(query, values)
        conn.commit()

    dur = (time.time()-t0)/60.0
    print(f'=======INGESTING gene_orthologs TOOK {dur:.2e} minutes=======')
