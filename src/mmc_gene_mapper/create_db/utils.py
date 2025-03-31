"""
Utilities for creating the gene mapper database
"""

import pathlib
import sqlite3


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
    create_gene_index(cursor)
    create_ensembl_index(cursor)
    create_ortholog_index(cursor)


def create_gene_index(cursor):
    idx_name = "symbol_idx"
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )
    cursor.execute(
        f"""
        CREATE INDEX {idx_name} ON NCBI_genes
        (species_taxon, symbol)
        """
    )

def create_ensembl_index(cursor):
    idx_name = "ncbi_idx"
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )
    cursor.execute(
        f"""
        CREATE INDEX {idx_name} ON NCBI_to_ENSEMBL
        (species_taxon, NCBI_id)
        """
    )

    idx_name = "ensembl_idx"
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )
    cursor.execute(
        f"""
        CREATE INDEX {idx_name} on NCBI_to_ENSEMBL
        (species_taxon, ENSEMBL_id)
        """
    )


def create_ortholog_index(cursor):
    idx_name = "ortholog_idx"
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )
    cursor.execute(
        f"""
        CREATE INDEX {idx_name} on orthologs
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
