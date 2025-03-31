import pathlib
import sqlite3


def ncbi_to_ensembl(
        db_path,
        ncbi_id_list,
        species_taxon,
        chunk_size=100,
        citation="NCBI"):
    does_path_exist(db_path)
    results = dict()
    with sqlite3.connect(db_path) as conn:

        citation = get_citation(
            conn=conn,
            name=citation
        )

        cursor = conn.cursor()
        for i0 in range(0, len(ncbi_id_list), chunk_size):
            values = [
                int(ii) for ii in ncbi_id_list[i0:i0+chunk_size]
            ]
            query=f"""
            SELECT
                NCBI_id,
                ENSEMBL_id
            FROM NCBI_to_ENSEMBL
            WHERE
                citation=?
            AND
                species_taxon=?
            AND
                NCBI_id IN {tuple(values)}
            """
            chunk = cursor.execute(
                query,
                (citation['idx'], species_taxon)).fetchall()
            for row in chunk:
                if row[0] not in results:
                    results[row[0]] = []
                results[row[0]].append(row[1])
    return results


def ensembl_to_ncbi(
        db_path,
        ensembl_id_list,
        species_taxon,
        chunk_size=100,
        citation="NCBI"):
    does_path_exist(db_path)
    results = dict()
    with sqlite3.connect(db_path) as conn:

        citation = get_citation(
            conn=conn,
            name=citation
        )

        cursor = conn.cursor()
        for i0 in range(0, len(ensembl_id_list), chunk_size):
            values = ensembl_id_list[i0:i0+chunk_size]
            query=f"""
            SELECT
                NCBI_id,
                ENSEMBL_id
            FROM NCBI_to_ENSEMBL
            WHERE
                citation=?
            AND
                species_taxon=?
            AND
                ENSEMBL_id IN {tuple(values)}
            """
            chunk = cursor.execute(
                query,
                (citation['idx'], species_taxon)).fetchall()
            for row in chunk:
                if row[0] not in results:
                    results[row[1]] = []
                results[row[1]].append(row[0])
    return results


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



def does_path_exist(db_path):
    db_path = pathlib.Path(db_path)
    if not db_path.is_file():
        raise RuntimeError(
            f"{db_path} is not a file"
        )
