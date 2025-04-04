import json
import pathlib
import sqlite3

import mmc_gene_mapper.create_db.metadata_tables as metadata_utils


def get_species_taxon(
        db_path,
        species_name,
        strict=False):
    does_path_exist(db_path)
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        result = _get_species_taxon(
            cursor=cursor,
            species_name=species_name,
            strict=strict)
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


def get_citation_from_bibliography(
        cursor,
        authority_idx,
        species_taxon):
    raw = cursor.execute(
        """
        SELECT
            citation
        FROM
            species_bibliography
        WHERE
            species_taxon=?
        AND
            authority=?
        """,
        (species_taxon, authority_idx)
    ).fetchall()

    results = set([r[0] for r in raw])

    if len(results) != 1:
        full_authority = cursor.execute(
            """
            SELECT
                name
            FROM authority
            WHERE id=?
            """,
            (authority,)
        )
        raise RuntimeError(
            f"There are {results} citations associated "
            f"with authority={full_authority}, "
            f"species_taxon={species_taxon}; "
            "unclear how to proceed"
        )
    citation_idx = results.pop()
    raw = cursor.execute(
        """
        SELECT
            name,
            metadata,
            id
        FROM citation
        WHERE id=?
        """,
        (citation_idx,)
    ).fetchall()
    assert len(raw) == 1
    return {
        "name": raw[0][0],
        "metadata": json.loads(raw[0][1]),
        "idx": raw[0][2]
    }


def get_authority_and_citation(
        conn,
        species_taxon,
        authority_name):
    """
    Return the full dicts characterizing authority
    and citaiton for a given (authority, species)
    combination

    Parameters
    ----------
    conn:
        sqlite3 connection
    species_taxon:
        an int identifying the species
    authority_name:
        a str; the name of the authority entry

    Returns
    -------
    A dict
        {'authority': full_authority,
         'citaiton': full_citation}

    Each value is a dict characterizing all the data
    carried around for that entity.

    Notes
    -----
    This function assumes there is one and only one
    valid citation for each (species_taxon, authority)
    combination. If more than one valid citation are found,
    an error will be thrown.
    """

    full_authority = metadata_utils.get_authority(
        conn=conn,
        name=authority_name
    )

    full_citation = get_citation_from_bibliography(
        cursor=conn.cursor(),
        authority_idx=full_authority['idx'],
        species_taxon=species_taxon
    )

    return {
        'authority': full_authority,
        'citation': full_citation
    }


def symbols_to_gene_identifiers(
        db_path,
        symbol_list,
        authority_name,
        species_taxon,
        chunk_size=100):

    results = {
        symbol: []
        for symbol in symbol_list
    }

    with sqlite3.connect(db_path) as conn:

        meta_source = get_authority_and_citation(
            conn=conn,
            authority_name=authority_name,
            species_taxon=species_taxon
        )

        full_authority = meta_source['authority']
        full_citation = meta_source['citation']

        authority_idx = full_authority["idx"]
        citation_idx = full_citation["idx"]

        cursor = conn.cursor()
        n_symbols = len(symbol_list)
        for i0 in range(0, n_symbols, chunk_size):
            values = symbol_list[i0:i0+chunk_size]
            n_values = len(values)

            query = """
                SELECT
                    symbol,
                    identifier
                FROM gene
                WHERE
                    citation=?
                AND
                    authority=?
                AND
                    species_taxon=?
                AND
                    symbol IN (
                """
            query += ",".join(['?']*n_values)
            query += ")"
            raw = cursor.execute(
                query,
                (citation_idx,
                 authority_idx,
                 species_taxon,
                 *values)
            )
            for row in raw:
                symbol = row[0]
                identifier = row[1]
                results[symbol].append(identifier)
        return {
            'metadata': {
                'authority': full_authority,
                'citation': full_citation
            },
            'mapping': results
        }


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
