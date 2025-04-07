"""
Utilities for creating the gene mapper database
"""

import json
import pathlib
import sqlite3
import time


def delete_index(cursor, idx_name):
    print(f'=======DELETING INDEX {idx_name}======')
    cursor.execute(
        f"""
        DROP INDEX IF EXISTS {idx_name}
        """
    )


def create_index(
        cursor,
        idx_name,
        table_name,
        column_tuple):

    delete_index(
        cursor=cursor,
        idx_name=idx_name)

    print(f'=======CREATING INDEX {idx_name}=======')

    cursor.execute(
        f"""
        CREATE INDEX {idx_name}
        ON {table_name}
        {tuple(column_tuple)}
        """
    )


def check_existence(db_path):
    db_path = pathlib.Path(db_path)
    if db_path.exists() and not db_path.is_file():
        raise RuntimeError(
            f"{db_path} exists but is not a file"
        )
    return db_path.exists()


data_table_list = (
    "gene",
    "gene_ortholog",
    "gene_equivalence",
)


metadata_table_list = (
    "authority",
    "citation"
)
