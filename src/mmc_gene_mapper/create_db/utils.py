"""
Utilities for creating the gene mapper database
"""

import pandas as pd
import pathlib
import sqlite3

import mmc_gene_mapper



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
        clobber=False):
    assert ensembl_path.is_file()
    assert ortholog_path.is_file()
    assert gene_info_path.is_file()

    db_path = pathlib.Path(db_path)
    if check_existence(db_path):
        if clobber:
            db_path.unlink()
        else:
            print(f'{db_path} exists; returning')
            return

    with sqlite3.connect(db_path) as conn:
        create_tables(conn)
        ingest_gene_info(
            conn=conn,
            data_path=gene_info_path
        )
        ingest_gene_to_ensembl(
            conn=conn,
            data_path=ensembl_path
        )
        ingest_orthologs(
            conn=conn,
            data_path=ortholog_path
        )
        create_indexes(conn)
            

def ingest_gene_info(conn, data_path):
    print('=======INGESTING GENE INFO=======')
    cursor = conn.cursor()
    chunk_size = 5000000
    query = """
    INSERT INTO NCBI_genes (
        species_taxon,
        NCBI_id,
        symbol
    ) VALUES (?, ?, ?)
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
                (int(row[0]), int(row[1]), row[2])
                for row in chunk
                if len(row) > 0
            ]
            if len(values) == 0:
                break
            print(f'    chunk {i0:.2e}, {i0+len(values):.2e}')
            i0 += len(chunk)
            cursor.executemany(query, values)
            conn.commit()


def ingest_gene_to_ensembl(conn, data_path):
    print('=======INGESTING GENE TO ENSEMBL=======')
    cursor = conn.cursor()
    data = pd.read_csv(data_path, delimiter='\t')
    n_rows = len(data)
    chunk_size = 50000
    query = """
    INSERT INTO NCBI_to_ENSEMBL (
        species_taxon,
        NCBI_id,
        ENSEMBL_id
    )
    VALUES (?, ?, ?)
    """
    for i0 in range(0, n_rows, chunk_size):
        chunk = data.iloc[i0:i0+chunk_size].to_dict(orient='records')
        values = [
            (int(row['#tax_id']),
             int(row['GeneID']),
             row['Ensembl_gene_identifier'])
            for row in chunk
        ]
        cursor.executemany(query, values)
        conn.commit()

def ingest_orthologs(conn, data_path):
    print('=======INGESTING ORTHOLOGS=======')
    cursor = conn.cursor()
    data = pd.read_csv(data_path, delimiter='\t')
    data = data[data['relationship'] == 'Ortholog']
    n_rows = len(data)
    chunk_size = 50000
    query = """
    INSERT INTO orthologs (
        species0,
        gene0,
        species1,
        gene1
    ) VALUES (?, ?, ?, ?)
    """
    for i0 in range(0, n_rows, chunk_size):
        chunk = data.iloc[i0:i0+chunk_size].to_dict(orient='records')
        values = [
            (int(row['#tax_id']),
             int(row['GeneID']),
             int(row['Other_tax_id']),
             int(row['Other_GeneID']))
            for row in chunk
        ]
        values += [
            (int(row['Other_tax_id']),
             int(row['Other_GeneID']),
             int(row['#tax_id']),
             int(row['GeneID']))
            for row in chunk
        ]
        cursor.executemany(query, values)
        conn.commit()


def create_tables(conn):
    print('=======CREATING TABLES=======')
    cursor = conn.cursor()

    cursor.execute(
        """
        CREATE TABLE NCBI_genes (
            species_taxon INTEGER,
            NCBI_id INTEGER,
            symbol STRING
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE NCBI_to_ENSEMBL (
            species_taxon INTEGER,
            NCBI_id INTEGER,
            ENSEMBL_id STRING
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE orthologs(
            species0 INTEGER,
            gene0 INTEGER,
            species1 INTEGER,
            gene1 INTEGER
        )
        """
    )

def create_indexes(conn):
    print('=======CREATING INDEXES=======')
    cursor = conn.cursor()
    conn.execute(
        """
        CREATE INDEX symbol_idx ON NCBI_genes
        (species_taxon, symbol)
        """
    )
    conn.execute(
        """
        CREATE INDEX ncbi_idx ON NCBI_to_ENSEMBL
        (species_taxon, NCBI_id)
        """
    )
    conn.execute(
        """
        CREATE INDEX ensembl_idx on NCBI_to_ENSEMBL
        (species_taxon, ENSEMBL_id)
        """
    )
    conn.execute(
        """
        CREATE INDEX ortholog_idx on orthologs
        (species0, species1, gene0)
        """
    )

def check_existence(db_path):
    db_path = pathlib.Path(db_path)
    if db_path.exists() and not db_path.is_file():
        raise RuntimeError(
            f"{db_path} exists but is not a file"
        )
    return db_path.exists()
