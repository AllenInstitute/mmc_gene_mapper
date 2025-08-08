import numpy as np
import pathlib
import sqlite3

import mmc_gene_mapper.utils.log_class as log_class
import mmc_gene_mapper.metadata.classes as metadata_classes
import mmc_gene_mapper.mapper.gene_characterization as gene_characterization


def detect_species_and_authority(
        db_path,
        gene_list,
        chunk_size=1000,
        guess_taxon=None,
        log=None):
    """
    Find the species and authority for a list of gene
    identifiers. Genes can be from an inhomogeneous list of
    authorities.

    Parameters
    ----------
    db_path:
        path to database being quried.
    gene_list:
        list of gene identifiers we are trying to map.
    chunk_size:
        an int; the number of genes to match to a species/authority
        at once
    guess_taxon:
       an optional int. This only comes into play if we end up
       having to choose a species based on gene symbols. If this
       is not None and there is a tie among taxons that match
       the provided gene symbols, this taxon will be chosen
       if it is one of the taxons participating in the tie.
    log:
        a logger class that implements an info()
        function (probably the CommandLog from cell_type_mapper)

    Returns
    --------
    A dict
        {
          "authority": np.array of strings indicating authority for genes
          "species": an instance of Species describing the species
        }

    Notes
    -----
    This function will first try to match species from NCBI and ENSEMBL
    identifiers. The function will assume that all genes in gene_list are
    from the same species. As such, it will accept the first match it finds
    as the truth. If it happens to find more than one species in the first
    chunk of data it analyzes, it will raise an exception.

    If ENSEMBL and NCBI identifiers are present, the function will
    try to match species to both. If one gives None and the other gives
    a non-None species, the non-None species will be preferred. If both
    give non-None answers and the answers differ, an exception will
    be thrown.

    If no species can be inferred from NCBI and ENSEMBL identifiers (or if
    no NCBI or ENSEMBL identifiers are present), the function will attempt
    to identify species using gene symbols by scanning all gene symbols and
    selecting the species that corresponds to the plurality of them.

    If no species is matched, authority will be set to a 'symbol' for all genes
    """

    if log is None:
        log = log_class.StdoutLog()

    gene_list = np.array(gene_list)

    # Do string matching to determine which genes
    # symbols, which are ENSEMBL IDs and which are
    # NCBI IDs

    authority = gene_characterization.characterize_gene_identifiers(
        db_path=db_path,
        gene_list=gene_list
    )
    authority = np.array(authority)

    symbol_idx = np.where(authority == 'symbol')[0]
    ensembl_idx = np.where(authority == 'ENSEMBL')[0]
    ncbi_idx = np.where(authority == 'NCBI')[0]

    n_found = (
        len(symbol_idx)
        + len(ensembl_idx)
        + len(ncbi_idx)
    )

    if n_found != len(gene_list):
        raise RuntimeError(
            "Not all genes accounted for in authorities "
            "('symbol', 'ENSEMBL', 'NCBI'); "
            f"input n: {len(gene_list)}\n"
            f"n authorities: {n_found}"
        )

    if len(ensembl_idx) > 0:
        ensembl_authority = _detect_species_and_authority(
            db_path=db_path,
            gene_list=gene_list[ensembl_idx],
            chunk_size=chunk_size
        )
    else:
        ensembl_authority = None

    if len(ncbi_idx) > 0:
        ncbi_authority = _detect_species_and_authority(
            db_path=db_path,
            gene_list=gene_list[ncbi_idx],
            chunk_size=chunk_size
        )
    else:
        ncbi_authority = None

    if ensembl_authority is None and ncbi_authority is None:
        species = None
    else:
        if ncbi_authority is None:
            species_name = ensembl_authority['species']
            species_taxon = ensembl_authority['species_taxon']
        elif ensembl_authority is None:
            species_name = ncbi_authority['species']
            species_taxon = ncbi_authority['species_taxon']
        else:
            ens = ensembl_authority['species_taxon']
            ncbi = ncbi_authority['species_taxon']
            if ens is None and ncbi is not None:
                species_name = ncbi_authority['species']
                species_taxon = ncbi_authority['species_taxon']
            elif ncbi is None and ens is not None:
                species_name = ensembl_authority['species']
                species_taxon = ensembl_authority['species_taxon']
            elif ens == ncbi:
                species_name = ncbi_authority['species']
                species_taxon = ncbi_authority['species_taxon']
            else:
                msg = (
                    "\nENSEMBL genes gave species "
                    f"'{ensembl_authority['species']}'"
                    "\nNCBI genes gave species "
                    f"'{ncbi_authority['species']}'"
                )
                raise InconsistentSpeciesError(msg)

        if species_name is None and species_taxon is None:
            species = None
        else:
            species = metadata_classes.Species(
                name=species_name,
                taxon=species_taxon
            )

    # if no species was identified, artificially set
    # authority == 'symbol'
    if species is None:
        authority = np.array(
            ['symbol']*len(authority)
        )

        # try to match a species by assuming all genes
        # are gene symbols
        symbol_species = _detect_species_from_symbols(
            gene_list=gene_list,
            db_path=db_path,
            chunk_size=chunk_size,
            guess_taxon=guess_taxon,
            log=log
        )

        if symbol_species is not None:
            species = metadata_classes.Species(
                name=symbol_species['species'],
                taxon=int(symbol_species['species_taxon'])
            )

    return {
        "authority": authority,
        "species": species
    }


def _detect_species_and_authority(
        db_path,
        gene_list,
        chunk_size=100):
    """
    Take a list of gene identifiers. Determine
    the authority and species associated
    with the genes.

    Parameters
    ----------
    db_path:
        path to the gene mapper database file to query
    gene_list:
        list of strings; the gene identifiers whose
        authority and species are being determined
    chunk_size:
        The number of gene identifiers to query at a time
        in a search for a match. Once a match is found,
        the corresponding (species, authority) pair will
        be assumed to be correct for all of the genes in
        gene_list. If an inconsistent answer is found in
        a single chunk, an exception will be raised.

    Returns
    -------
    A dict
        {
          "authority": "name of the authority associated "
                       "with the genes (ENSEMBL or NCBI)"

          "species": "the name of the species associated "
                     "with the genes"

          "species_taxon": the integer ID of the species
        }

    Notes
    -----
    This function can be very slow as it has to query the
    database on an unindexed string (the gene identifier)
    to find a species. If no species is found, this could
    involve many costly queries to the database. Try not
    to run it on lists of gene symbols (which will not
    be found in the database).

    If the genes do not exist in the gene table, a "None"
    will be returned for all fields in the output.
    """

    bad_result = {
        "authority": None,
        "species": None,
        "species_taxon": None
    }

    n_genes = len(gene_list)
    with sqlite3.connect(db_path) as conn:
        for i0 in range(0, n_genes, chunk_size):
            sub_list = gene_list[i0:i0+chunk_size]
            query = """
                SELECT
                    gene.species_taxon,
                    authority.name
                FROM gene
                JOIN authority on authority.id = gene.authority
                WHERE
                    gene.identifier IN (
            """
            query += ",".join(["?"]*len(sub_list))
            query += ")"
            cursor = conn.cursor()
            results = cursor.execute(
                query,
                sub_list
            ).fetchall()

            if len(results) == 0:
                continue

            authority_set = sorted(set(r[1] for r in results))
            species_set = sorted(set(r[0] for r in results))
            error_msg = ""
            if len(authority_set) > 1:
                error_msg += (
                    f"\nMultiple authorities inferred: {authority_set}"
                )
            if len(species_set) > 1:
                error_msg += f"\nMultiple species inferred : {species_set}"
            if len(error_msg) > 0:
                error_msg = (
                    "Could not infer species and authority from gene list."
                    f"{error_msg}"
                )
                raise InconsistentSpeciesError(error_msg)

            species_query = """
            SELECT
                name
            FROM NCBI_species
            WHERE
                id=?
            LIMIT 1
            """
            species_name = cursor.execute(
                species_query,
                species_set
            ).fetchall()

            return {
                "authority": authority_set[0],
                "species": species_name[0][0],
                "species_taxon": species_set[0]
            }

    return bad_result


def _detect_species_from_symbols(
        db_path,
        gene_list,
        chunk_size,
        guess_taxon=None,
        log=None):
    """
    Parameters
    ----------
    db_path:
        path to the gene mapper database file to query
    gene_list:
        list of strings; the gene symbols whose
        species are being determined
    chunk_size:
        The number of gene identifiers to query at a time
        in a search for a match. Because of the ambiguity
        in gene symbols, this function scans through all
        of the provided symbols to make sure a consistent
        species can be found.
    guess_taxon:
       an optional int. If not None and the symbols match
       multiple species equally well, choose this one.
    log:
        a logger class that implements an info()
        function (probably the CommandLog from cell_type_mapper)

    Returns:
    --------
    A dict
        {"species": the name of the matching species
         "species_taxon": the taxon ID of the matching species}

    Notes
    -----
    Most gene symbols will match to several species. This function
    will query all of the provided symbols, finding all species that
    have genes corresponding to those symbols. Whichever species
    occurs most frequently will be the chosen species.

    If no species match, this function will return None.

    If more than one species match, this function will raise an
    InconsistentSpeciesError

    If guess_taxon is not None and another species is a better
    match for the gene symbols, return that better match, anyway.
    """
    if log is None:
        log = log_class.StdoutLog()

    n_genes = len(gene_list)
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        raw_lookup = dict()
        for i0 in range(0, n_genes, chunk_size):
            chunk = gene_list[i0: i0+chunk_size]
            query = """
            SELECT
                symbol,
                species_taxon
            FROM
                gene
            WHERE symbol in (
            """
            query += ",".join(["?"]*len(chunk))
            query += ")"
            raw_results = cursor.execute(query, chunk).fetchall()

            # first create a mapping from gene symbol to
            # a set of unique taxons, to guard against the
            # same (symbol, taxon) pair being ingested more
            # than once (perhaps, from multiple authorities or
            # citations)
            for row in raw_results:
                symbol = row[0]
                taxon = row[1]
                if symbol not in raw_lookup:
                    raw_lookup[symbol] = set()
                raw_lookup[symbol].add(taxon)

        if len(raw_lookup) == 0:
            return None

        # Now tally the votes
        taxon_candidates = np.concat(
            [sorted(raw_lookup[gene])
             for gene in sorted(raw_lookup.keys())]
        )
        taxon_candidates, vote_counts = np.unique(
            taxon_candidates,
            return_counts=True
        )

        max_dex = np.argmax(vote_counts)
        max_votes = vote_counts[max_dex]
        chosen_max = (vote_counts == max_votes)
        chosen_taxon = taxon_candidates[chosen_max]

        if len(chosen_taxon) > 1:
            if guess_taxon is not None:
                if guess_taxon in chosen_taxon:
                    msg = (
                        "Input gene symbols are consistent "
                        f"with several species taxons. Taxon {guess_taxon} "
                        "is one of them. Will assume the genes are "
                        "aligned to that taxon."
                    )
                    log.warn(msg)
                    chosen_taxon = [guess_taxon]

        chosen_taxon_name = []
        for taxon_id in chosen_taxon:
            species_name_query = """
            SELECT
                name
            FROM NCBI_species
            WHERE
                id=?
            LIMIT 1
            """
            name_result = cursor.execute(
                species_name_query,
                (int(taxon_id),)
            ).fetchall()
            chosen_taxon_name.append(name_result[0][0])

    if len(chosen_taxon) > 1:
        msg = (
            "The gene symbols you gave are consistent with "
            "more than one species: "
        )
        msg += ", ".join(
            [f"{name}:{taxon}"
             for name, taxon in zip(chosen_taxon_name, chosen_taxon)]
        )
        raise InconsistentSpeciesError(msg)

    return {
        "species": chosen_taxon_name[0],
        "species_taxon": chosen_taxon[0]
    }


def detect_if_genes(
        db_path,
        gene_list,
        chunk_size=1000):
    """
    Iterate over a list of genes. Return True if at least one
    of the genes is a valid gene identifier or symbol. Return
    False otherwise.
    """
    if db_path is None:
        raise ValueError(
            "you passed db_path = None to "
            "mmc_gene_mapper.mapper.species_dection.detect_if_genes; "
            "must specify a db_path"
        )
    db_path = pathlib.Path(db_path)
    if not db_path.is_file():
        raise ValueError(
            f"{db_path} is not a file "
            "(error in "
            "mmc_gene_mapper.mapper.species_detection.detect_if_genes)"
        )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        for i0 in range(0, len(gene_list), chunk_size):
            chunk = gene_list[i0:i0+chunk_size]
            id_query = """
            SELECT
                COUNT(identifier)
            FROM
                gene
            WHERE
                identifier IN (
            """
            id_query += ",".join([f'"{g}"' for g in chunk])
            id_query += ")"
            result = cursor.execute(id_query).fetchall()
            if result[0][0] > 0:
                return True
            symbol_query = """
            SELECT
                COUNT(symbol)
            FROM
                gene
            WHERE
                symbol IN (
            """
            symbol_query += ",".join([f'"{g}"' for g in chunk])
            symbol_query += ")"
            result = cursor.execute(symbol_query).fetchall()
            if result[0][0] > 0:
                return True
    return False


class InconsistentSpeciesError(Exception):
    pass
