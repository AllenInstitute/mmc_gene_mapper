"""
Define utilities around metadata tables
(authority and citation)
"""
import json
import pathlib
import sqlite3
import time

import mmc_gene_mapper.create_db.utils as db_utils


def create_metadata_tables(conn):
    cursor = conn.cursor()

    cursor.execute(
        """
        CREATE TABLE citation(
            id INTEGER,
            name STR,
            metadata STRING
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE authority(
           id INTEGER,
           name STR
        )
        """
    )


def get_citation(
        conn,
        name,
        strict=False):
    """
    Query the database at conn for a specific citation.

    Parameters
    ----------
    conn:
        the sqlite3 connection to the database
    name:
        a str; the name of the citation being requested
    strict:
        if True and there is no match, raise a ValueError.
        if False and there is no match, return None

    Returns
    -------
    A dict
        {
            "name": the name of the citation,
            "metadata": the metadata associated with the citation
            "idx": the integer index of the citation in the database
        }
    """
    cursor = conn.cursor()
    return _get_citation(
        cursor=cursor,
        name=name,
        strict=strict
    )


def _get_citation(
        cursor,
        name,
        strict=False):
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
        if strict:
            raise ValueError(
                f"No citation corresponding to name {name}"
            )
        return None

    return {
        "name": results[0][0],
        "idx": results[0][1],
        "metadata": json.loads(results[0][2])
    }



def delete_citation(conn, name):
    _delete_metadata(
        conn=conn,
        table_name="citation",
        name=name
    )


def insert_citation(
        conn,
        name,
        metadata_dict):

    pre_existing = get_citation(
        conn=conn,
        name=name
    )

    if pre_existing is not None:
        raise ValueError(
            f"citation name {name} already exists"
        )

    cursor = conn.cursor()

    metadata_str = json.dumps(metadata_dict)

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
            raise ValueError(
                f"citation {name} already exists; "
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


def get_authority(
        conn,
        name,
        strict=False):

    cursor = conn.cursor()
    return _get_authority(
        cursor=cursor,
        name=name,
        strict=strict)


def _get_authority(
        cursor,
        name,
        strict=False):

    results = cursor.execute(
        """
        SELECT name, id
        FROM authority
        WHERE name=?
        """,
        (name,)
    ).fetchall()

    if len(results) > 1:
        raise ValueError(
            f"More than one authority corresponding to {name}"
        )

    if len(results) == 0:
        if strict:
            raise ValueError(
                f"No authority corresponding to name {name}"
            )
        return None

    return {
        "name": results[0][0],
        "idx": results[0][1]
    }



def delete_authority(conn, name):
    _delete_metadata(
        conn=conn,
        table_name="authority",
        name=name
    )


def insert_authority(
        conn,
        name):

    pre_existing = get_authority(
        conn=conn,
        name=name)

    if pre_existing is not None:
        return pre_existing

    cursor = conn.cursor()

    auth_max = cursor.execute(
        """
        SELECT MAX(id) FROM authority
        """
    ).fetchall()[0][0]

    if auth_max is None:
        auth_idx = 0
    else:
        auth_idx = auth_max + 1

    cursor.execute(
        """
        INSERT INTO authority(
            id,
            name
        )
        VALUES (?, ?)
        """,
        (auth_idx, name)
    )
    return auth_idx


def insert_unique_authority(
        conn,
        name,
        clobber=False):
    """
    If citation already exists, delete it.
    """
    pre_existing_citation = get_authority(
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
            delete_authority(
                conn=conn,
                name=name
            )

    auth_idx = insert_authority(
        conn=conn,
        name=name
    )
    return auth_idx


def _delete_metadata(conn, table_name, name):
    t0 = time.time()
    print(f'=======DELETING {table_name} {name}=======')
    cursor = conn.cursor()

    pre_existing = cursor.execute(
        f"""
        SELECT
            id
        FROM
            {table_name}
        WHERE name=?
        """,
        (name,)
    ).fetchall()

    if len(pre_existing) != 1:
        raise ValueError(
            f"{len(pre_existing)} {table_name} with name {name}; "
            "expected exactly one"
        )

    target_idx = pre_existing[0][0]

    for data_table_name in db_utils.data_table_list:
        cursor.execute(
            f"""
            DELETE FROM {data_table_name}
            WHERE citation=?
            """,
            (target_idx,)
        )
        print(f"    DELETED FROM {table_name}")

    cursor.execute(
        f"""
        DELETE FROM {table_name}
        WHERE id=?
        """,
        (target_idx, )
    )

    conn.commit()
    dur = time.time() - t0
    print(f"    DELETING TOOK {dur:.2e} seconds")

