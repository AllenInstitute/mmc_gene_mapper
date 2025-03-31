"""
Utilities for creating the gene mapper database
"""

import json
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
            symbol STRING,
            citation INTEGER
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE NCBI_to_ENSEMBL (
            species_taxon INTEGER,
            NCBI_id INTEGER,
            ENSEMBL_id STRING,
            citation INTEGER
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE NCBI_orthologs(
            species0 INTEGER,
            gene0 INTEGER,
            species1 INTEGER,
            gene1 INTEGER,
            citation INTEGER
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE citation(
            id INTEGER,
            name STR,
            metadata STRING
        )
        """
    )

def create_indexes(conn):
    print('=======CREATING INDEXES=======')
    cursor = conn.cursor()
    create_gene_index(cursor)
    create_ensembl_index(cursor)
    create_ncbi_ortholog_index(cursor)


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
        (citation, species_taxon, symbol)
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
        (citation, species_taxon, NCBI_id)
        """
    )

    idx_name = "ncbi_ensembl_idx"
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )
    cursor.execute(
        f"""
        CREATE INDEX {idx_name} on NCBI_to_ENSEMBL
        (citation, species_taxon, ENSEMBL_id)
        """
    )


def create_ncbi_ortholog_index(cursor):
    idx_name = "ncbi_ortholog_idx"
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )
    cursor.execute(
        f"""
        CREATE INDEX {idx_name} on NCBI_orthologs
        (citation, species0, species1, gene0)
        """
    )


def insert_citation(
        conn,
        name,
        metadata_dict):

    cursor = conn.cursor()

    metadata_str = json.dumps(metadata_dict)

    pre_existing = cursor.execute(
        """
        SELECT COUNT(*) FROM citation
        WHERE name = ?
        """,
        (name,)
    ).fetchall()

    if pre_existing[0][0] != 0:
        raise ValueError(
            f"citation name {name} already exists"
        )

    citation_max = cursor.execute(
        """
        SELECT MAX(id) FROM citation
        """
    ).fetchall()[0][0]

    if citation_max is None:
        citation_idx = 0
    else:
        citation_idx = citation_max + 1

    cursor.execute(
        """
        INSERT INTO citation(
            id,
            name,
            metadata
        )
        VALUES (?, ?, ?)
        """,
        (citation_idx, name, metadata_str)
    )
    return citation_idx


def check_existence(db_path):
    db_path = pathlib.Path(db_path)
    if db_path.exists() and not db_path.is_file():
        raise RuntimeError(
            f"{db_path} exists but is not a file"
        )
    return db_path.exists()
