import numpy as np
import pathlib
import sqlite3

import mmc_gene_mapper.utils.log_class as log_class
import mmc_gene_mapper.utils.str_utils as str_utils
import mmc_gene_mapper.query_db.query as query_utils


def detect_species_and_authority(
        db_path,
        gene_list,
        chunk_size=1000,
        guess_taxon=None,
        log=None,
        clean_ensembl=True):
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
    clean_ensembl:
        if True, loop over gene_list to remove version suffixes
        from ENSEMBL IDs

    Returns
    --------
    A dict
        {
          "authority": np.array of strings indicating authority for genes
          "species": an instance of Species describing the species
        }

    Notes
    -----
    This function will query the gene table, attempting to match the
    genes in gene_list to known gene identifiers and gene symbols.

    It will do a poll of the results, assuming your data is aligned to
    the species with the most matches (if multiple species tie, an error
    will be raised).

    It will then loop back through the genes, assigning them to the
    authority that matches the chosen species.

    Any genes that are not assigned an authority will be assumed to be
    gene symbols. This means that ENSEMBL IDs and NCBI IDs that have not
    been ingested into the database will be treated as gene symbols.
    """

    if log is None:
        log = log_class.StdoutLog()

    if clean_ensembl:
        gene_list = str_utils.remove_ensembl_versions(gene_list)

    gene_list = np.array(gene_list)
    n_genes = len(gene_list)

    from_identifier = species_from_identifier(
        db_path=db_path,
        gene_list=gene_list,
        log=log
    )

    from_symbol = species_from_symbol(
        db_path=db_path,
        gene_list=gene_list,
        log=log
    )

    votes = []

    for gene in gene_list:
        these_votes = []
        if gene in from_identifier:
            these_votes += [
                val['species_taxon']
                for val in from_identifier[gene]
            ]
        if gene in from_symbol:
            these_votes += [
                val['species_taxon']
                for val in from_symbol[gene]
            ]
        if len(these_votes) == 0:
            continue
        these_votes, these_cts = np.unique(these_votes, return_counts=True)
        max_ct = these_cts.max()
        valid = np.where(these_cts == max_ct)
        votes.append(these_votes[valid])

    if len(votes) == 0:
        chosen_species = None
    else:
        votes = np.concatenate(votes)
        votes, cts = np.unique(votes, return_counts=True)
        chosen_cts = cts.max()
        chosen_dex = np.where(cts == chosen_cts)
        chosen_taxon = votes[chosen_dex]

        chosen_species = None
        with sqlite3.connect(db_path) as conn:
            cursor = conn.cursor()
            species_lookup = {
                ii: query_utils.get_species(
                    cursor=cursor,
                    species=int(ii)
                )
                for ii in chosen_taxon
            }

            if len(chosen_taxon) > 1:

                # try to break the species degeneracy with the
                # guess_taxon
                broke_degeneracy = False
                msg = (
                    f"{chosen_cts} of your genes were consistent with "
                    "the following species:\n"
                )
                for val in species_lookup.values():
                    msg += f"'{val}'\n"

                if guess_taxon is not None:
                    if guess_taxon in chosen_taxon:
                        broke_degeneracy = True
                        chosen_species = species_lookup[guess_taxon]
                        msg += (
                            "using guess to resolve degeneracy "
                            "in favor of '{speces_lookup[guess_taxon]}'"
                        )
                        log.warn(msg)

                if not broke_degeneracy:
                    msg += "Unable to break this degeneracy"
                    raise InconsistentSpeciesError(msg)

            if chosen_species is None:
                chosen_species = query_utils.get_species(
                    cursor=cursor,
                    species=int(chosen_taxon[0])
                )

        log.info(
            f"Based on {chosen_cts} genes, your input data is from species "
            f"'{chosen_species}'"
        )

    # assign authorities to the genes
    authority = np.array([None]*n_genes)

    if chosen_species is not None:
        for i_gene, gene in enumerate(gene_list):
            assigned = False
            if gene in from_identifier:
                for val in from_identifier[gene]:
                    if val['species_taxon'] == chosen_species.taxon:
                        authority[i_gene] = val['authority']
                        assigned = True
                        break
            if not assigned:
                if gene in from_symbol:
                    authority[i_gene] = 'symbol'

    # any authorities that are still unknown will be assumed to be
    # symbols
    unknown_idx = np.array(
        [ii for ii in range(n_genes) if authority[ii] is None]
    )
    if len(unknown_idx) > 0:
        authority[unknown_idx] = 'symbol'

    return {
        "authority": authority,
        "species": chosen_species
    }


def species_from_identifier(
        db_path,
        gene_list,
        chunk_size=10000,
        log=None):
    """
    Return a dict mapping the gene_identifiers in gene_list
    to (authority, species_taxon) pairs
    """
    return _species_from_column(
        db_path=db_path,
        gene_list=gene_list,
        col_name="identifier",
        chunk_size=chunk_size,
        log=log
    )


def species_from_symbol(
        db_path,
        gene_list,
        chunk_size=10000,
        log=None):
    """
    Return a dict mapping the gene_identifiers in gene_list
    to (authority, species_taxon) pairs
    """
    return _species_from_column(
        db_path=db_path,
        gene_list=gene_list,
        col_name="symbol",
        chunk_size=chunk_size,
        log=log
    )


def _species_from_column(
        db_path,
        gene_list,
        col_name,
        chunk_size,
        log=None):

    if log is None:
        log = log_class.StdoutLog()

    if col_name not in ("identifier", "symbol"):
        raise RuntimeError(
            f"{col_name} is not a valid col_name; "
            "must be either 'identifier' or 'symbol'"
        )

    db_path = pathlib.Path(db_path)
    if not db_path.is_file():
        raise RuntimeError(f"db_path {db_path} is not a file")

    n_genes = len(gene_list)
    mapping = dict()
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        for i0 in range(0, n_genes, chunk_size):
            chunk = gene_list[i0: i0+chunk_size]
            query = f"""
            SELECT
                gene.{col_name},
                gene.species_taxon,
                authority.name
            FROM
                gene
            JOIN authority on gene.authority = authority.id
            WHERE
                gene.{col_name} IN (
            """
            query += ",".join(["?"]*len(chunk))
            query += ")"

            # to prevent multiple citation entries from skewing the count
            query += f"""
                GROUP BY gene.{col_name}, gene.authority, gene.species_taxon
            """

            results = cursor.execute(
                query,
                chunk
            ).fetchall()
            for row in results:
                if row[0] not in mapping:
                    mapping[row[0]] = []
                mapping[row[0]].append(
                    {'authority': row[2],
                     'species_taxon': row[1]}
                )
    return mapping


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


class MultipleAuthorityError(Exception):
    pass
