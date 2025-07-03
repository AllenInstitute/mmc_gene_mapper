import json
import pathlib
import sqlite3

import mmc_gene_mapper.metadata.classes as metadata_classes
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils


def get_species(
        cursor,
        species):
    """
    Return a metadata_classes.Species

    Parameters
    ----------
    cursor:
        a sqlite3 cursor
    species:
        either a string or an int specifying the species

    Returns
    -------
    a metadata_classes.Species containing the taxon and
    name of the specified species
    """
    if isinstance(species, int):
        species_taxon = species
        species_name = _get_species_name(
            cursor=cursor,
            species_taxon=species_taxon
        )
    elif isinstance(species, str):
        species_name = species
        species_taxon = _get_species_taxon(
            cursor=cursor,
            species_name=species_name
        )
    else:
        raise ValueError(
            f"Cannot infer species from '{species}'"
            f"of type {type(species)}"
        )

    return metadata_classes.Species(
        taxon=species_taxon,
        name=species_name
    )


def _get_species_taxon(
        cursor,
        species_name):
    """
    Return the integer identifier for a species

    Parameters
    ----------
    cursor:
        a sqlite3.cursor
    species_name:
        human-readable name of the species

    Returns
    -------
    an integer; the taxon ID of the species
    """
    try:
        species_taxon = int(species_name)
        return species_taxon
    except ValueError:
        pass

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
        raise ValueError(
            f"{len(results)} species match name {species_name}\n"
            f"{results}"
        )
    elif len(results) == 0:
        raise ValueError(
            f"no species match for {species_name}"
        )
    return results[0][0]


def _get_species_name(
        cursor,
        species_taxon):
    """
    Return the name of a species

    Parameters
    ----------
    cursor:
        a sqlite3.cursor
    species_taxon:
        an int; the species taxon

    Returns
    -------
    a string; the name of the species
    """
    results = cursor.execute(
       """
       SELECT
           name
       FROM
           NCBI_species
       WHERE
           id=?
       """,
       (species_taxon,)
    ).fetchall()

    if len(results) == 0:
        raise ValueError(
            f"no species match for {species_taxon}"
        )
    return results[0][0]


def get_citation_from_bibliography(
        cursor,
        authority_idx,
        species,
        require_symbols=False):
    """
    Return the full entry for a citation based on
    authority and species (assumes there is only one
    valid citation).

    Parameters
    ----------
    cursor:
        sqlite3 cursor object
    authority_idx:
        the integer index of the relevant authority
    species:
        an instance of Species
    require_symbols:
        if True, only consider citation, authority pairs
        with symbols associated with them

    Returns
    -------
    A dict representing the citation that associates
    that species with that authority (via the gene table)
        {
            "name": name of the citation,
            "idx": numerical index of the citation,
            "metadata": dict containing citation's metadata
        }

    Notes
    -----
    If require_symbols is False and multiple citations are
    returned, try running again with require_symbols=True
    to see if you only get one citation back.
    """
    query = """
        SELECT
            citation
        FROM
            species_bibliography
        WHERE
            species_taxon=?
        AND
            authority=?
    """
    if require_symbols:
        query += """
            AND
               has_symbols=1
            """

    raw = cursor.execute(
        query,
        (species.taxon, authority_idx)
    ).fetchall()

    results = set([r[0] for r in raw])

    if len(results) != 1:
        if not require_symbols:
            return get_citation_from_bibliography(
                cursor=cursor,
                authority_idx=authority_idx,
                species=species,
                require_symbols=True)
        else:
            full_authority = cursor.execute(
                """
                SELECT
                    name
                FROM authority
                WHERE id=?
                """,
                (authority_idx,)
            ).fetchall()

            if len(full_authority) == 0:
                full_authority = authority_idx

            raise ValueError(
                f"There are {len(results)} citations associated "
                f"with authority={full_authority}, "
                f"species_taxon={species.taxon}; "
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
        species,
        authority_name,
        require_symbols=False):
    """
    Return the full dicts characterizing authority
    and citation for a given (authority, species)
    combination

    Parameters
    ----------
    conn:
        sqlite3 connection
    species:
        an instance of Species
    authority_name:
        a str; the name of the authority entry
    require_symbols:
        if True, only consider (citation, authority) pairs
        with associated gene symbols

    Returns
    -------
    A dict
        {'authority': full_authority,
         'citation': full_citation}

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

    if full_authority is None:
        raise ValueError(
            f"authority {authority_name} does not exist "
            "in this database"
        )

    full_citation = get_citation_from_bibliography(
        cursor=conn.cursor(),
        authority_idx=full_authority['idx'],
        species=species,
        require_symbols=require_symbols
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
        species,
        chunk_size=100):

    if src_column not in ('symbol', 'id', 'identifier'):
        raise RuntimeError(
            f"{src_column} not a valid src column for gene table"
        )

    if dst_column not in ('symbol', 'id', 'identifier'):
        raise RuntimeError(
            f"{dst_column} not a valid dst column for gene table"
        )

    if src_column == 'symbol':
        require_symbols = True
    else:
        require_symbols = False

    results = {
        val: set()
        for val in src_list
    }

    if src_column == 'symbol':
        mapping_src = metadata_classes.Authority('symbol')
    else:
        mapping_src = metadata_classes.Authority(authority_name)

    if dst_column == 'symbol':
        mapping_dst = metadata_classes.Authority('symbol')
    else:
        mapping_dst = metadata_classes.Authority(authority_name)

    with sqlite3.connect(db_path) as conn:

        meta_source = get_authority_and_citation(
            conn=conn,
            authority_name=authority_name,
            species=species,
            require_symbols=require_symbols
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
                 species.taxon,
                 *values)
            ).fetchall()
            for row in raw:
                val = row[0]
                identifier = row[1]
                results[val].add(identifier)

        for key in results:
            results[key] = sorted(results[key])

        return {
            'metadata': metadata_classes.MappingMetadata(
                src=mapping_src,
                dst=mapping_dst,
                citation=full_citation
            ).serialize(),
            'mapping': results
        }


def get_equivalent_genes_from_identifiers(
        db_path,
        input_authority_name,
        output_authority_name,
        input_gene_list,
        species,
        citation_name,
        chunk_size=100):

    does_path_exist(db_path)

    id_translation = translate_gene_identifiers(
        db_path=db_path,
        src_column='identifier',
        dst_column='id',
        src_list=input_gene_list,
        authority_name=input_authority_name,
        species=species,
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
        species_taxon=species.taxon,
        citation_name=citation_name,
        chunk_size=chunk_size
    )

    mapping_dict = mapping_dict_to_identifiers(
        db_path=db_path,
        mapping_dict=equivalence['mapping'],
        key_authority_name=input_authority_name,
        key_species=species,
        value_authority_name=output_authority_name,
        value_species=species
    )

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
            query = """
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
        'metadata': metadata_classes.MappingMetadata(
            src=metadata_classes.Authority(input_authority_name),
            dst=metadata_classes.Authority(output_authority_name),
            citation=full_citation
        ).serialize(),
        'mapping': results
    }


def get_ortholog_genes_from_identifiers(
        db_path,
        authority_name,
        src_species,
        dst_species,
        src_gene_list,
        citation_name,
        chunk_size=100):

    does_path_exist(db_path)

    id_translation = translate_gene_identifiers(
        db_path=db_path,
        src_column='identifier',
        dst_column='id',
        src_list=src_gene_list,
        authority_name=authority_name,
        species=src_species,
        chunk_size=chunk_size
    )

    id_values = set()
    for key in id_translation['mapping']:
        id_values = id_values.union(set(id_translation['mapping'][key]))

    orthologs = _get_ortholog_genes(
        db_path=db_path,
        authority_name=authority_name,
        src_species_taxon=src_species.taxon,
        src_genes=sorted(id_values),
        dst_species_taxon=dst_species.taxon,
        citation_name=citation_name,
        chunk_size=chunk_size
    )

    mapping_dict = mapping_dict_to_identifiers(
        db_path=db_path,
        mapping_dict=orthologs['mapping'],
        key_authority_name=authority_name,
        key_species=src_species,
        value_authority_name=authority_name,
        value_species=dst_species
    )

    unmapped_genes = set(src_gene_list)-set(mapping_dict.keys())
    for gene in unmapped_genes:
        mapping_dict[gene] = []

    return {
        'metadata': orthologs['metadata'],
        'mapping': mapping_dict
    }


def _get_ortholog_genes(
        db_path,
        authority_name,
        src_species_taxon,
        src_genes,
        dst_species_taxon,
        citation_name,
        chunk_size=50):
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
    src_genes:
        list of ints indicating the NCBI IDs of the
        specified genes
    dst_species_taxon:
        int indicating the species in which you
        want to find ortholog genes
    citation_name:
        string indicating the source of the ortholog
        assignments you want to use
    chunk_size:
        how many values to process at once
    """
    does_path_exist(db_path)

    results = dict()

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
            gene_chunk = tuple(src_genes[i0:i0+chunk_size])
            query = """
                SELECT
                    gene,
                    ortholog_group
                FROM
                    gene_ortholog
                WHERE
                    authority=?
                AND
                    citation=?
                AND
                    species=?
                AND
                    gene IN (
            """
            query += ",".join(['?']*len(gene_chunk))
            query += ")"
            gene_to_ortholog = cursor.execute(
                query,
                (authority_idx,
                 citation['idx'],
                 src_species_taxon,
                 *gene_chunk)
            ).fetchall()
            n_raw = len(gene_to_ortholog)

            gene_to_ortholog = {
                row[0]: row[1]
                for row in gene_to_ortholog
            }
            if len(gene_to_ortholog) != n_raw:
                raise ValueError(
                    "gene_to_ortholog was not unique"
                )

            ortholog_chunk = tuple(
                set(gene_to_ortholog.values())
            )

            query = """
                SELECT
                    ortholog_group,
                    gene
                FROM
                    gene_ortholog
                WHERE
                    authority=?
                AND
                    citation=?
                AND
                    species=?
                AND
                    ortholog_group IN (
            """
            query += ",".join(['?']*len(ortholog_chunk))
            query += ")"
            ortholog_to_other_gene = cursor.execute(
                query,
                (authority_idx,
                 citation['idx'],
                 dst_species_taxon,
                 *ortholog_chunk)
            ).fetchall()
            n_raw = len(ortholog_to_other_gene)

            ortholog_to_other_gene = {
                row[0]: row[1]
                for row in ortholog_to_other_gene
            }
            if len(ortholog_to_other_gene) != n_raw:
                raise ValueError(
                    "ortholog_to_other_gene was not unique"
                )

            for gene in gene_to_ortholog:
                orth = gene_to_ortholog[gene]
                if orth not in ortholog_to_other_gene:
                    continue
                if gene not in results:
                    results[gene] = []
                results[gene].append(ortholog_to_other_gene[orth])

    return {
        'metadata': {
            'citation': citation['metadata']
        },
        'mapping': results
    }


def mapping_dict_to_identifiers(
        db_path,
        mapping_dict,
        key_authority_name,
        key_species,
        value_authority_name,
        value_species):
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
    key_species:
        Species associated with the keys
    value_authority_name:
        name of the authority defining the values
    value_species:
        Species associated with the values

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
        species=key_species
    )
    error_msg += key_error
    if len(error_msg) > 0:
        raise MappingError(error_msg)

    value_list = set()
    for key in mapping_dict:
        value_list = value_list.union(set(mapping_dict[key]))
    value_list = sorted(value_list)

    value_mapping, _ = _strict_mapping_from_id(
        db_path=db_path,
        value_list=list(value_list),
        authority_name=value_authority_name,
        species=value_species
    )

    new_dict = dict()
    for key in mapping_dict:
        new_key = key_mapping[key][0]
        new_dict[new_key] = []
        for val in mapping_dict[key]:
            if val in value_mapping:
                if len(value_mapping[val]) == 1:
                    new_dict[new_key].append(
                        value_mapping[val][0]
                    )
                elif len(value_mapping[val]) > 1:
                    raise RuntimeError(
                        f"More than one mapping value for {new_key}"
                    )

    return new_dict


def _strict_mapping_from_id(
        db_path,
        value_list,
        authority_name,
        species):

    return _strict_mapping(
        db_path=db_path,
        src_column='id',
        dst_column='identifier',
        src_list=value_list,
        authority_name=authority_name,
        species=species,
        allow_none=False
    )


def _strict_mapping(
        db_path,
        src_column,
        dst_column,
        src_list,
        authority_name,
        species,
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
        species=species
    )['mapping']

    error_msg = ""
    for val in src_list:
        if len(mapping[val]) != 1:
            if len(mapping[val]) == 0 and allow_none:
                pass
            else:
                error_msg += (
                   f"id: {val} authority: {authority_name} "
                   f"species: {species.taxon} "
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
