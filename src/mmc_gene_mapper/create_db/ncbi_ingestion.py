import json
import pandas as pd
import pathlib
import sqlite3

import mmc_gene_mapper
import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.create_db.utils as db_utils


def create_database(db_name, clobber=False):

    data_dir = pathlib.Path(mmc_gene_mapper.__file__).parent / 'data'
    db_dir = pathlib.Path(mmc_gene_mapper.__file__).parent / 'db_files'
    assert data_dir.is_dir()
    assert db_dir.is_dir()
    gene_info_path = data_dir / 'gene_info'
    ortholog_path = data_dir / 'gene_orthologs'
    ensembl_path = data_dir / 'gene2ensembl'

    db_path = db_dir / db_name

    _create_database(
        db_path=db_path,
        gene_info_path=gene_info_path,
        ensembl_path=ensembl_path,
        ortholog_path=ortholog_path,
        clobber=clobber
    )

def _create_database(
        db_path,
        gene_info_path,
        ortholog_path,
        ensembl_path,
        clobber=False,
        citation_tag='NCBI'):
    assert ensembl_path.is_file()
    assert ortholog_path.is_file()
    assert gene_info_path.is_file()

    db_path = pathlib.Path(db_path)
    if db_utils.check_existence(db_path):
        if clobber:
            db_path.unlink()
        else:
            print(f'{db_path} exists; returning')
            return

    metadata_dict = {
        'downloaded_on': timestamp.get_timestamp()
    }

    with sqlite3.connect(db_path) as conn:
        db_utils.create_tables(conn)

        citation_idx = db_utils.insert_citation(
            conn=conn,
            name='NCBI',
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
    with open(data_path, 'r') as src:
        header = src.readline()
        while True:
            chunk = [
                src.readline().strip().split()[:3]
                for ii in range(chunk_size)
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

