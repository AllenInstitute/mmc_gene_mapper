import json
import numpy as np
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


def translate_gene_identifiers(
        db_path,
        src_column,
        dst_column,
        src_list,
        authority_name,
        species_taxon,
        chunk_size=100):

    if src_column not in ('symbol', 'id', 'identifier'):
        raise RuntimeError(
            f"{src_column} not a valid src column for gene table"
        )

    if dst_column not in ('symbol', 'id', 'identifier'):
        raise RuntimeError(
            f"{dst_column} not a valid dst column for gene table"
        )

    results = {
        val: []
        for val in src_list
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
        n_vals = len(src_list)
        for i0 in range(0, n_vals, chunk_size):
            values = src_list[i0:i0+chunk_size]
            n_values = len(values)

            query = f"""
                SELECT
                    {src_column},
                    {dst_column}
                FROM gene
                WHERE
                    citation=?
                AND
                    authority=?
                AND
                    species_taxon=?
                AND
                    {src_column} IN (
                """
            query += ",".join(['?']*n_values)
            query += ")"
            raw = cursor.execute(
                query,
                (citation_idx,
                 authority_idx,
                 species_taxon,
                 *values)
            ).fetchall()
            for row in raw:
                val = row[0]
                identifier = row[1]
                results[val].append(identifier)
        return {
            'metadata': {
                'authority': full_authority,
                'citation': full_citation
            },
            'mapping': results
        }


def get_equivalent_genes_from_identifiers(
        db_path,
        input_authority_name,
        output_authority_name,
        input_gene_list,
        species_name,
        citation_name,
        chunk_size=100):

    species_taxon = get_species_taxon(
        db_path=db_path,
        species_name=species_name,
        strict=True
    )

    id_translation = translate_gene_identifiers(
        db_path=db_path,
        src_column='identifier',
        dst_column='id',
        src_list=input_gene_list,
        authority_name=input_authority_name,
        species_taxon=species_taxon,
        chunk_size=chunk_size
    )

    id_values = set()
    for key in id_translation['mapping']:
        id_values = id_values.union(set(id_translation['mapping'][key]))

    equivalence = _get_equivalent_genes(
        db_path=db_path,
        input_id_list=sorted(id_values),
        input_authority_name=input_authority_name,
        output_authority_name=output_authority_name,
        species_taxon=species_taxon,
        citation_name=citation_name,
        chunk_size=chunk_size
    )
    
    mapping_dict = mapping_dict_to_identifiers(
        db_path=db_path,
        mapping_dict=equivalence['mapping'],
        key_authority_name=input_authority_name,
        key_species_taxon=species_taxon,
        value_authority_name=output_authority_name,
        value_species_taxon=species_taxon
    )

    equivalence['metadata']['key_authority'] = input_authority_name
    equivalence['metadata']['value_authority'] = output_authority_name

    unmapped_genes = set(input_gene_list)-set(mapping_dict.keys())
    for gene in unmapped_genes:
        mapping_dict[gene] = []

    return {
        'metadata': equivalence['metadata'],
        'mapping': mapping_dict
    }


def _get_equivalent_genes(
        db_path,
        input_id_list,
        input_authority_name,
        output_authority_name,
        species_taxon,
        chunk_size=100,
        citation_name="NCBI"):
    """
    This takes as input ids (the integer IDs used in the database)
    """

    input_id_list = [int(ii) for ii in input_id_list]

    does_path_exist(db_path)
    results = {
        ii: []
        for ii in input_id_list
    }
    with sqlite3.connect(db_path) as conn:

        full_citation = metadata_utils.get_citation(
            conn=conn,
            name=citation_name
        )

        input_auth = metadata_utils.get_authority(
            conn=conn,
            name=input_authority_name
        )["idx"]

        output_auth = metadata_utils.get_authority(
            conn=conn,
            name=output_authority_name
        )["idx"]

        cursor = conn.cursor()
        for i0 in range(0, len(input_id_list), chunk_size):
            values = [
                ii for ii in input_id_list[i0:i0+chunk_size]
            ]
            n_values = len(values)
            query="""
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
                gene0 IN (
            """
            query += ",".join(['?']*n_values)
            query += ")"
            chunk = cursor.execute(
                query,
                (full_citation['idx'],
                 input_auth,
                 output_auth,
                 species_taxon,
                 *values)).fetchall()
            for row in chunk:
                results[row[0]].append(row[1])

    return {
        'metadata': {
            'key_authority': input_auth,
            'value_authority': output_auth,
            'citation': full_citation
        },
        'mapping': results
    }


def get_orthologs(
        db_path,
        authority_name,
        src_species_taxon,
        src_genes,
        dst_species_taxon,
        citation_name):
    """
    Parameters
    ----------
    db_path:
        path to the database file
    authority_name:
        string indicating according to what authority
        (NCBI or ENSEMBL) we want orthologs
    src_species_taxon:
        int indicating the species of the
        specified genes
    dst_genes:
        list of ints indicating the NCBI IDs of the
        sepcified genes
    dst_species_taxon:
        int indicating the species in which you
        want to find ortholog genes
    citation_name:
        string indicating the source of the ortholog
        assignments you want to use
    """
    does_path_exist(db_path)

    results = dict()
    chunk_size = 500

    with sqlite3.connect(db_path) as conn:
        citation = metadata_utils.get_citation(
            conn=conn,
            name=citation_name
        )

        authority_idx = metadata_utils.get_authority(
            conn=conn,
            name=authority_name
        )["idx"]

        cursor = conn.cursor()

        for i0 in range(0, len(src_genes), chunk_size):
            values = tuple(src_genes[i0:i0+chunk_size])
            query = """
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
                    gene0 IN (
                """
            query += ",".join(['?']*len(values))
            query += ")"
            raw = cursor.execute(
                query,
                (authority_idx,
                 citation['idx'],
                 src_species_taxon,
                 dst_species_taxon,
                 *values)
            )
            for row in raw:
                if row[0] not in results:
                    results[row[0]] = []
                results[row[0]].append(row[1])

    return {
        'metadata': {
            'citation': json.loads(citation['metadata']),
            'src_species_taxon': src_species_taxon,
            'dst_species_taxon': dst_species_taxon
        },
        'mapping': results
    }


def mapping_dict_to_identifiers(
        db_path,
        mapping_dict,
        key_authority_name,
        key_species_taxon,
        value_authority_name,
        value_species_taxon):
    """
    Take a dict mapping gene ids (ints) to other
    gene ids. Convert it to a dict mapping
    gene identifiers (strings) to other gene identifiers.

    Parameters
    ----------
    db_path:
        path to the database
    mapping_dict:
        keys should be ints; values should be lists of ints
    key_authority_name:
        name of the authority defining the keys
    key_species_taxon:
        int defining the species associated with the keys
    value_authority_name:
        name of the authority defining the values
    value_species_taxon:
        int defining the species associated with the values

    Returns
    -------
    A dict
       equivalent to mapping dict, but id integers
       have been converted to identifier strings
    """

    error_msg = ""
    key_mapping, key_error = _strict_mapping_from_id(
        db_path=db_path,
        value_list=list(mapping_dict.keys()),
        authority_name=key_authority_name,
        species_taxon=key_species_taxon
    )
    error_msg += key_error

    value_list = set()
    for key in mapping_dict:
        value_list = value_list.union(set(mapping_dict[key]))
    value_list = sorted(value_list)

    value_mapping, value_error = _strict_mapping_from_id(
        db_path=db_path,
        value_list=list(value_list),
        authority_name=value_authority_name,
        species_taxon=value_species_taxon
    )
    error_msg += value_error
    if len(error_msg) > 0:
        raise MappingError(error_msg)

    new_dict = dict()
    for key in mapping_dict:
        new_key = key_mapping[key][0]
        new_dict[new_key] = []
        for val in mapping_dict[key]:
            new_dict[new_key].append(
                value_mapping[val][0]
            )

    return new_dict


def _strict_mapping_from_id(
        db_path,
        value_list,
        authority_name,
        species_taxon):

    return _strict_mapping(
        db_path=db_path,
        src_column='id',
        dst_column='identifier',
        src_list=value_list,
        authority_name=authority_name,
        species_taxon=species_taxon,
        allow_none=False
    )


def _strict_mapping(
        db_path,
        src_column,
        dst_column,
        src_list,
        authority_name,
        species_taxon,
        allow_none=True):
    """
    demand 1:1 mapping
    """

    mapping = translate_gene_identifiers(
        db_path=db_path,
        src_column='id',
        dst_column='identifier',
        src_list=src_list,
        authority_name=authority_name,
        species_taxon=species_taxon
    )['mapping']

    error_msg = ""
    for val in src_list:
        if len(mapping[val]) != 1:
            if len(mapping[val]) == 0 and allow_none:
                pass
            else:
                error_msg += (
                   f"id: {val} authority: {authority_name} "
                    f"species: {species_taxon} "
                    f"n: {len(mapping[val])}\n"
                )
    return mapping, error_msg


def does_path_exist(db_path):
    db_path = pathlib.Path(db_path)
    if not db_path.is_file():
        raise RuntimeError(
            f"{db_path} is not a file"
        )


class MappingError(Exception):
    pass
