import gzip
import json
import pandas as pd
import pathlib
import sqlite3
import time

import mmc_gene_mapper
import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.create_db.utils as db_utils
import mmc_gene_mapper.query_db.query as db_query
import mmc_gene_mapper.download.ftp_utils as ftp_utils


def update_ncbi_data(
        db_name,
        data_dir):

    t0 = time.time()
    data_dir = pathlib.Path(data_dir)
    if data_dir.exists():
        if not data_dir.is_dir():
            raise RuntimeError(
                f"{data_dir} is not a dir"
            )
    else:
        data_dir.mkdir(parents=True)


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
        db_name=db_name,
        data_dir=data_dir,
        clobber=True
    )
    dur = (time.time()-t0)/60.0
    print(f'SUCCESS; WHOLE PROCESS TOOK {dur:.2e} minutes')


def ingest_ncbi_data(
        db_name,
        data_dir,
        clobber=False):

    data_dir = pathlib.Path(data_dir)
    db_dir = pathlib.Path(mmc_gene_mapper.__file__).parent / 'db_files'
    assert data_dir.is_dir()
    assert db_dir.is_dir()
    gene_info_path = data_dir / 'gene_info.gz'
    ortholog_path = data_dir / 'gene_orthologs.gz'
    ensembl_path = data_dir / 'gene2ensembl.gz'
    metadata_path = data_dir / 'metadata.json'

    db_path = db_dir / db_name

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
            db_utils.create_tables(conn)

        pre_existing_citation = db_query.get_citation(
           conn=conn,
           name=citation_name
        )

        if pre_existing_citation is not None:
            if not clobber:
                raise RuntimeError(
                    f"citation {citation_name} already exists; "
                    "run with clobber=True to overwrite"
                )
            else:
                db_utils.delete_citation(
                    conn=conn,
                    name=citation_name
                )

        citation_idx = db_utils.insert_citation(
            conn=conn,
            name=citation_name,
            metadata_dict=metadata_dict
        )

        ingest_gene_info(
            conn=conn,
            data_path=gene_info_path,
            citation_idx=citation_idx
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
        db_utils.create_indexes(conn)
            

def ingest_gene_info(conn, data_path, citation_idx):
    print('=======INGESTING GENE INFO=======')
    data_path = pathlib.Path(data_path)
    cursor = conn.cursor()
    chunk_size = 5000000
    query = """
    INSERT INTO NCBI_genes (
        species_taxon,
        NCBI_id,
        symbol,
        citation
    ) VALUES (?, ?, ?, ?)
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
                (int(row[0]), int(row[1]), row[2], citation_idx)
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
    cursor = conn.cursor()
    data = pd.read_csv(data_path, delimiter='\t')
    n_rows = len(data)
    chunk_size = 50000
    query = """
    INSERT INTO NCBI_to_ENSEMBL (
        species_taxon,
        NCBI_id,
        ENSEMBL_id,
        citation
    )
    VALUES (?, ?, ?, ?)
    """
    for i0 in range(0, n_rows, chunk_size):
        chunk = data.iloc[i0:i0+chunk_size].to_dict(orient='records')
        values = [
            (int(row['#tax_id']),
             int(row['GeneID']),
             row['Ensembl_gene_identifier'],
             citation_idx)
            for row in chunk
        ]
        cursor.executemany(query, values)
        conn.commit()


def ingest_orthologs(
        conn,
        data_path,
        citation_idx):
    print('=======INGESTING ORTHOLOGS=======')

    cursor = conn.cursor()
    data = pd.read_csv(data_path, delimiter='\t')
    data = data[data['relationship'] == 'Ortholog']
    n_rows = len(data)
    chunk_size = 50000
    query = """
    INSERT INTO NCBI_orthologs (
        species0,
        gene0,
        species1,
        gene1,
        citation
    ) VALUES (?, ?, ?, ?, ?)
    """
    for i0 in range(0, n_rows, chunk_size):
        chunk = data.iloc[i0:i0+chunk_size].to_dict(orient='records')
        values = [
            (int(row['#tax_id']),
             int(row['GeneID']),
             int(row['Other_tax_id']),
             int(row['Other_GeneID']),
             citation_idx)
            for row in chunk
        ]
        values += [
            (int(row['Other_tax_id']),
             int(row['Other_GeneID']),
             int(row['#tax_id']),
             int(row['GeneID']),
             citation_idx)
            for row in chunk
        ]
        cursor.executemany(query, values)
        conn.commit()
