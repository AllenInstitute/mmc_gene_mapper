import json
import pathlib
import sqlite3


def get_species_taxon(
        db_path,
        species_name):
    does_path_exist(db_path)
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        results = cursor.execute(
            """
            SELECT
                id
            FROM
                NCBI_species
            WHERE
                name=?
            """,
            (species_name,)
        ).fetchall()
    if len(results) > 1:
        raise RuntimeError(
            f"{len(results)} species match name {name}\n"
            f"{results}"
        )
    elif len(results) == 0:
        return None
    return results[0][0]


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


def get_ncbi_orthologs(
        db_path,
        src_species,
        src_genes,
        dst_species,
        citation):
    """
    Parameters
    ----------
    db_path:
        path to the database file
    src_species:
        string or int indicating the species of the
        specified genes
    dst_genes:
        list of ints indicating the NCBI IDs of the
        sepcified genes
    dst_species:
        string or int indicating the species in which you
        want to find ortholog genes
    citation:
        string indicating the source of the ortholog
        assignments you want to use
    """
    does_path_exist(db_path)

    src_taxon = species_to_taxon(
        db_path=db_path,
        species=src_species)

    dst_taxon = species_to_taxon(
        db_path=db_path,
        species=dst_species)

    results = dict()
    chunk_size = 500

    with sqlite3.connect(db_path) as conn:
        citation = get_citation(
            conn=conn,
            name=citation
        )

        cursor = conn.cursor()

        for i0 in range(0, len(src_genes), chunk_size):
            values = tuple(src_genes[i0:i0+chunk_size])
            raw = cursor.execute(
                f"""
                SELECT
                    gene0,
                    gene1
                FROM
                    NCBI_orthologs
                WHERE
                    citation=?
                AND
                    species0=?
                AND
                    species1=?
                AND
                    gene0 IN {values}
                """,
                (citation['idx'], src_taxon, dst_taxon)
            )
            for row in raw:
                if row[0] not in results:
                    results[row[0]] = []
                results[row[0]].append(row[1])

    return {
        'metadata': {
            'provenance': json.loads(citation['metadata']),
            'src_species_taxon': src_taxon,
            'dst_species_taxon': dst_taxon
        },
        'mapping': results
    }


def species_to_taxon(db_path, species):
    if isinstance(species, int):
        return species

    try:
        result = int(species)
    except ValueError:
        result = get_species_taxon(
            db_path=db_path,
            species_name=species
        )
    return result


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
