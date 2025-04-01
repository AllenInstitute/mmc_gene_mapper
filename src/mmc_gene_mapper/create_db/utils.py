"""
Utilities for creating the gene mapper database
"""

import json
import pathlib
import sqlite3
import time


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

    idx_name = "gene_species_idx"
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )
    cursor.execute(
        f"""
        CREATE INDEX {idx_name} ON NCBI_genes
        (citation, species_taxon, NCBI_id)
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


def insert_unique_citation(
        conn,
        name,
        metadata_dict,
        clobber=False):
    """
    If citation already exists, delete it.
    """
    pre_existing_citation = get_citation(
       conn=conn,
       name=name
    )

    if pre_existing_citation is not None:
        if not clobber:
            raise RuntimeError(
                f"citation {citation_name} already exists; "
                "run with clobber=True to overwrite"
            )
        else:
            delete_citation(
                conn=conn,
                name=name
            )

    citation_idx = insert_citation(
        conn=conn,
        name=name,
        metadata_dict=metadata_dict
    )
    return citation_idx


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


def delete_citation(conn, name):
    t0 = time.time()
    print(f'=======DELETING CITATION {name}=======')
    cursor = conn.cursor()
    pre_existing = cursor.execute(
        """
        SELECT
            id
        FROM
            citation
        WHERE name=?
        """,
        (name,)
    ).fetchall()

    if len(pre_existing) != 1:
        raise ValueError(
            f"{len(pre_existing)} citations with name {name}; "
            "expected exactly one"
        )

    target_idx = pre_existing[0][0]

    for idx_name in ("symbol_idx",
                     "ncbi_idx",
                     "ncbi_ensembl_idx",
                     "ncbi_ortholog_idx"):
        cursor.execute(
            f"DROP INDEX IF EXISTS {idx_name}"
        )

    print("    DELETED INDEXES")

    for table_name in ('NCBI_genes',
                       'NCBI_to_ENSEMBL',
                       'NCBI_orthologs'):
        cursor.execute(
            f"""
            DELETE FROM {table_name}
            WHERE citation=?
            """,
            (target_idx,)
        )
        print(f"    DELETED FROM {table_name}")

    cursor.execute(
        """
        DELETE FROM citation
        WHERE id=?
        """,
        (target_idx, )
    )

    conn.commit()
    dur = time.time() - t0
    print(f"    DELETING TOOK {dur:.2e} seconds")


def get_citation(conn, name):

    cursor = conn.cursor()

    results = cursor.execute(
        """
        SELECT name, id, metadata
        FROM citation
        WHERE name=?
        """,
        (name,)
    ).fetchall()

    if len(results) > 1:
        raise ValueError(
            f"More than one citation corresponding to {name}"
        )

    if len(results) == 0:
        return None

    return {
        "name": results[0][0],
        "idx": results[0][1],
        "metadata": results[0][2]
    }


def check_existence(db_path):
    db_path = pathlib.Path(db_path)
    if db_path.exists() and not db_path.is_file():
        raise RuntimeError(
            f"{db_path} exists but is not a file"
        )
    return db_path.exists()
