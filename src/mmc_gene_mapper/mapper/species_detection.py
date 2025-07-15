import copy
import numpy as np
import sqlite3

import mmc_gene_mapper.utils.str_utils as str_utils
import mmc_gene_mapper.metadata.classes as metadata_classes


def detect_species_and_authority(
        db_path,
        gene_list,
        chunk_size=1000):
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

    Returns
    --------
    A dict
        {
          "authority": np.array of strings indicating authority for genes
          "species": an instance of Species describing the species
        }

    Notes
    -----
    This function assumes that all genes in gene_list are from the same
    species. As such, it will accept the first match it finds as the
    truth. If it happens to find more than one species in the first chunk
    of data it analyzes, it will raise an exception. If you pass it a list
    of gene symbols, "species" and "species_taxon" will be None
    """

    gene_list = np.array(gene_list)

    # Do string matching to determine which genes
    # symbols, which are ENSEMBL IDs and which are
    # NCBI IDs

    authority = str_utils.characterize_gene_identifiers(
        gene_list
    )
    authority = np.array(authority)

    symbol_idx = np.where(authority=='symbol')[0]
    ensembl_idx = np.where(authority=='ENSEMBL')[0]
    ncbi_idx = np.where(authority=='NCBI')[0]

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
            if ens != ncbi:
                msg = (
                    f"\nENSEMBL genes gave species '{ensembl_authority['species']}'"
                    f"\nNCBI genes gave species '{ncbi_authority['species']}'"
                )
                raise InconsistentSpeciesError(msg)
            species_name = ncbi_authority['species']
            species_taxon = ncbi_authority['species_taxon']

        species = metadata_classes.Species(
            name=species_name,
            taxon=species_taxon
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


class InconsistentSpeciesError(Exception):
    pass
