import pathlib
import sqlite3


def ncbi_to_ensembl(
        db_path,
        ncbi_id_list,
        species_taxon,
        chunk_size=100):
    does_path_exist(db_path)
    results = dict()
    with sqlite3.connect(db_path) as conn:

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
                species_taxon=?
            AND
                NCBI_id IN {tuple(values)}
            """
            chunk = cursor.execute(query, (species_taxon, )).fetchall()
            for row in chunk:
                if row[0] not in results:
                    results[row[0]] = []
                results[row[0]].append(row[1])
    return results
                


def does_path_exist(db_path):
    db_path = pathlib.Path(db_path)
    if not db_path.is_file():
        raise RuntimeError(
            f"{db_path} is not a file"
        )
