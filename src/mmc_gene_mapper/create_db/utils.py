"""
Utilities for creating the gene mapper database
"""
import pathlib

import mmc_gene_mapper.utils.file_utils as file_utils


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

    if len(column_tuple) > 1:
        cursor.execute(
            f"""
            CREATE INDEX {idx_name}
            ON {table_name}
            {tuple(column_tuple)}
            """
        )
    else:
        cursor.execute(
            f"""
            CREATE INDEX {idx_name}
            ON {table_name}
            ({column_tuple[0]})
            """
        )


def check_existence(db_path):
    """
    Check that the file at db_path exists.
    If it exists but is not a file, raise a
    NotAFileError

    (this function exists to prevent us from
    creating a sqlite3 connection to a new file when
    we do not intend to)
    """
    db_path = pathlib.Path(db_path)
    if db_path.exists() and not db_path.is_file():
        raise file_utils.NotAFileError(
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
