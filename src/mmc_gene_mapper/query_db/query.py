import json
import pathlib
import sqlite3


import mmc_gene_mapper.create_db.metadata_tables as metadata_utils


def get_species_taxon(
        db_path,
        species_name):
    does_path_exist(db_path)
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        result = _get_species_taxon(
            cursor=cursor,
            species_name=species_name)
    return result

def _get_species_taxon(
        cursor,
        species_name,
        strict=False):
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
        if strict:
            raise RuntimeError(
                f"no species match for {species_name}"
            )
        return None
    return results[0][0]


def get_equivalent_genes(
        db_path,
        input_id_list,
        input_authority,
        output_authority,
        species_taxon,
        chunk_size=100,
        citation="NCBI"):

    does_path_exist(db_path)
    results = dict()
    with sqlite3.connect(db_path) as conn:

        citation = metadata_utils.get_citation(
            conn=conn,
            name=citation
        )

        input_auth = metadata_utils.get_authority(
            conn=conn,
            name=input_authority
        )["idx"]

        output_auth = metadata_utils.get_authority(
            conn=conn,
            name=output_authority
        )["idx"]

        cursor = conn.cursor()
        for i0 in range(0, len(input_id_list), chunk_size):
            values = [
                int(ii) for ii in input_id_list[i0:i0+chunk_size]
            ]
            query=f"""
            SELECT
                gene0,
                gene1
            FROM gene_equivalence
            WHERE
                citation=?
            AND
                authority0=?
            AND
                authority1=?
            AND
                species_taxon=?
            AND
                gene0 IN {tuple(values)}
            """
            chunk = cursor.execute(
                query,
                (citation['idx'],
                 input_auth,
                 output_auth,
                 species_taxon)).fetchall()
            for row in chunk:
                if row[0] not in results:
                    results[row[0]] = []
                results[row[0]].append(row[1])
    return results


def get_orthologs(
        db_path,
        authority,
        src_species,
        src_genes,
        dst_species,
        citation):
    """
    Parameters
    ----------
    db_path:
        path to the database file
    authority:
        string indicating according to what authority
        (NCBI or ENSEMBL) we want orthologs
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
        citation = metadata_utils.get_citation(
            conn=conn,
            name=citation
        )

        authority_idx = metadata_utils.get_authority(
            conn=conn,
            name=authority
        )["idx"]

        cursor = conn.cursor()

        for i0 in range(0, len(src_genes), chunk_size):
            values = tuple(src_genes[i0:i0+chunk_size])
            raw = cursor.execute(
                f"""
                SELECT
                    gene0,
                    gene1
                FROM
                    gene_ortholog
                WHERE
                    authority=?
                AND
                    citation=?
                AND
                    species0=?
                AND
                    species1=?
                AND
                    gene0 IN {values}
                """,
                (authority_idx,
                 citation['idx'],
                 src_taxon,
                 dst_taxon)
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


def does_path_exist(db_path):
    db_path = pathlib.Path(db_path)
    if not db_path.is_file():
        raise RuntimeError(
            f"{db_path} is not a file"
        )
