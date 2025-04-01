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
import mmc_gene_mapper.download.ftp_utils as ftp_utils


def update_ncbi_data(
        db_path,
        data_dir,
        do_download=True):

    t0 = time.time()
    data_dir = pathlib.Path(data_dir)
    if data_dir.exists():
        if not data_dir.is_dir():
            raise RuntimeError(
                f"{data_dir} is not a dir"
            )
    else:
        data_dir.mkdir(parents=True)


    if do_download:
        host = 'ftp.ncbi.nlm.nih.gov'
        mapping = {
            'gene/DATA/gene_info.gz': data_dir/'gene_info.gz',
            'gene/DATA/gene2ensembl.gz': data_dir/'gene2ensembl.gz',
            'gene/DATA/gene_orthologs.gz': data_dir/'gene_orthologs.gz'
        }
        metadata_dst = data_dir/'metadata.json'

        ftp_utils.download_files_from_ftp(
            ftp_host=host,
            file_dst_mapping=mapping,
            metadata_dst=metadata_dst
        )

    ingest_ncbi_data(
        db_path=db_path,
        data_dir=data_dir,
        clobber=True
    )
    dur = (time.time()-t0)/60.0
    print(f'SUCCESS; WHOLE PROCESS TOOK {dur:.2e} minutes')


def ingest_ncbi_data(
        db_path,
        data_dir,
        clobber=False):

    data_dir = pathlib.Path(data_dir)
    assert data_dir.is_dir()

    gene_info_path = data_dir / 'gene_info.gz'
    ortholog_path = data_dir / 'gene_orthologs.gz'
    ensembl_path = data_dir / 'gene2ensembl.gz'
    metadata_path = data_dir / 'metadata.json'

    _ingest_ncbi_data(
        db_path=db_path,
        gene_info_path=gene_info_path,
        ensembl_path=ensembl_path,
        ortholog_path=ortholog_path,
        metadata_path=metadata_path,
        clobber=clobber,
        citation_name='NCBI'
    )

def _ingest_ncbi_data(
        db_path,
        gene_info_path,
        ortholog_path,
        ensembl_path,
        metadata_path,
        clobber=False,
        citation_name='NCBI'):

    file_utils.assert_is_file(ensembl_path)
    file_utils.assert_is_file(ortholog_path)
    file_utils.assert_is_file(gene_info_path)
    file_utils.assert_is_file(metadata_path)

    db_exists = False
    db_path = pathlib.Path(db_path)
    if db_utils.check_existence(db_path):
        db_exists = True

    with open(metadata_path, 'rb') as src:
        metadata_dict = json.load(src)

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
            authority_idx=auth_idx
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

        data_utils.create_data_indexes(conn)
            

def ingest_gene_info(conn, data_path, authority_idx):
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
        identifier
    ) VALUES (?, ?, ?, ?, ?)
    """
    i0 = 0
    if data_path.suffix == '.gz':
        open_fn = gzip.open
        mode = 'rb'
    else:
        open_fn = open
        mode = 'r'

    with open_fn(data_path, mode=mode) as src:
        header = src.readline()
        while True:
            chunk = [
                src.readline().strip().split()[:3]
                for ii in range(chunk_size)
            ]
            if mode == 'rb':
                chunk = [
                    [el.decode('utf-8') for el in row]
                    for row in chunk
                ]

            values = [
                (authority_idx,
                 int(row[0]),
                 int(row[1]),
                 row[2],
                 f'NCBIGene:{row[1]}')
                for row in chunk
                if len(row) > 0
            ]

            if len(values) == 0:
                break

            print(f'    chunk {i0:.2e}, {i0+len(values):.2e}')

            i0 += len(chunk)
            cursor.executemany(query, values)
            conn.commit()


def ingest_gene_to_ensembl(conn, data_path, citation_idx):
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
    chunk_size = 50000

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

    for i0 in range(0, n_rows, chunk_size):
        chunk = data.iloc[i0:i0+chunk_size].to_dict(orient='records')
        values = []
        for row in chunk:
            try:
                ensembl_gene_id = str_utils.int_from_identifier(
                    row['Ensembl_gene_identifier']
                )
            except ValueError:
                continue

            ncbi_gene_id = int(row['GeneID'])
            taxon_id = int(row['#tax_id'])

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


def ingest_orthologs(
        conn,
        data_path,
        citation_idx):
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
