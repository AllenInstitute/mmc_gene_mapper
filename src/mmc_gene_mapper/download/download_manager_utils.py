"""
Functions supporting the download manager
"""

import pathlib
import sqlite3

import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.utils.file_utils as file_utils


def create_download_db(db_path):
    """
    Create a database at db_path and initialize the downloads table.

    Raise an exception of db_path already exists.
    """

    db_path = pathlib.Path(db_path)
    if db_path.exists():
        raise ValueError(f"{db_path} already exists")

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE downloads (
                host STRING,
                src_path STRING,
                local_path STRING,
                hash STRING,
                downloaded_on STRING
            )
            """
        )


def remove_record(db_path, host, src_path):
    db_path = pathlib.Path(db_path)
    file_utils.assert_is_file(db_path)
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            DELETE FROM downloads
            WHERE
                host=?
            AND
                src_path=?
            """,
            (host, src_path)
        )


def insert_record(
        db_path,
        host,
        src_path,
        local_path):

    db_path = pathlib.Path(db_path)
    file_utils.assert_is_file(db_path)
    file_utils.assert_is_file(local_path)
    hash_val = file_utils.hash_from_path(local_path)
    date_val = timestamp.get_timestamp()
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.executemany(
            """
            INSERT INTO downloads (
                host,
                src_path,
                local_path,
                hash,
                downloaded_on
            )
            VALUES (?, ?, ?, ?, ?)
            """,
            [(host, src_path, local_path, hash_val, date_val)]
        )


def get_record(
        db_path,
        host,
        src_path):
    file_utils.assert_is_file(db_path)
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        results = cursor.execute(
            """
            SELECT
                host,
                src_path,
                local_path,
                hash,
                downloaded_on
            FROM downloads
            WHERE
                host=?
            AND
                src_path=?
            """,
            (host, src_path)
        )

    return [
        {'host': r[0],
         'src_path': r[1],
         'local_path': r[2],
         'hash': r[3],
         'downloaded_on': r[4]}
        for r in results
    ]
