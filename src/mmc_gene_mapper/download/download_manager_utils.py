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
    """
    Delete all the records in the downloads table corresponding
    to a specific (host, src) pair.

    Parameters
    ----------
    db_path:
        path to the database file being affected
    host:
        a string; the value in the 'host' field corresponding
        to the records to be deleted
    src_path:
        a string; the value in the src_path field corresponding
        to the records being deleted

    Returns
    -------
    None
        appropriate rows in the downloads table are deleted
    """
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
    """
    Insert a row representing a specific data file
    into the downloads table

    Parameters
    ----------
    db_path:
        path to the database file being updated
    host:
        a string; the host from which the data file
        was downloaded
    src_path:
        a string; the path on host from which the
        data file was downloaded
    local_path:
        the path on the local system to which the
        data file was downloaded

    Returns
    -------
    None
        database is updated with the appropriate
        information, including the current timestamp
        and the hash of the data file (which are
        calculated by this function)

    Notes
    -----
    if local_path is not a valid file, an exception is
    raised
    """

    db_path = pathlib.Path(db_path)
    file_utils.assert_is_file(db_path)
    file_utils.assert_is_file(local_path)
    hash_val = file_utils.hash_from_path(local_path)
    date_val = timestamp.get_timestamp()
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
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
            (host, src_path, local_path, hash_val, date_val)
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
