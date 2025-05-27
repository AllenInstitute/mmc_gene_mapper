"""
utility functions for gene mapper class
"""

import copy
import numpy as np
import sqlite3

import mmc_gene_mapper.utils.str_utils as str_utils
import mmc_gene_mapper.create_db.utils as db_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.ncbi_ingestion as ncbi_ingestion
import mmc_gene_mapper.create_db.species_ingestion as species_ingestion
import mmc_gene_mapper.create_db.bkbit_ingestion as bkbit_ingestion
import mmc_gene_mapper.create_db.ortholog_ingestion as ortholog_ingestion


def create_mapper_database(
        db_path,
        download_manager,
        tmp_dir,
        force_download,
        data_file_spec=None):

    with sqlite3.connect(db_path) as conn:
        data_utils.create_data_tables(conn)
        metadata_utils.create_metadata_tables(conn)

    ncbi_ingestion.ingest_ncbi_data(
        db_path=db_path,
        download_manager=download_manager,
        clobber=False
    )

    species_ingestion.ingest_species_data(
        db_path=db_path,
        download_manager=download_manager,
        tmp_dir=tmp_dir,
        force_download=force_download
    )

    if data_file_spec is not None:
        for file_spec in data_file_spec:
            if file_spec['type'] == 'bkbit':
                bkbit_ingestion.ingest_bkbit_genes(
                    db_path=db_path,
                    bkbit_path=file_spec['path']
                )

    if data_file_spec is not None:
        for file_spec in data_file_spec:
            if file_spec['type'] == 'bkbit':
                continue
            elif file_spec['type'] == 'hmba_orthologs':
                ortholog_ingestion.ingest_hmba_orthologs(
                    db_path=db_path,
                    hmba_file_path=file_spec['path'],
                    citation_name=file_spec['name'],
                    clobber=False
                )
            else:
                raise RuntimeError(
                    f"cannot parse file of type {file_spec['type']}"
                )

    # only after all data has been ingested
    with sqlite3.connect(db_path) as conn:
        data_utils.create_data_indexes(conn)


def create_bibliography_table(
        db_path):
    """
    Create species_bibliography which lists all
    combinations of authority, citation, species in the
    gene table.
    """
    print("=======CREATING BIBLIOGRAPHY TABLE=======")
    table_name = "species_bibliography"
    index_name = "species_bibliography_idx"
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        db_utils.delete_index(
            cursor=cursor,
            idx_name=index_name
        )
        cursor.execute(
            f"""
            DROP TABLE IF EXISTS {table_name}
            """
        )
        cursor.execute(
            f"""
            CREATE TABLE {table_name} (
                authority INTEGER,
                citation INTEGER,
                species_taxon INTEGER,
                has_symbols INTEGER
            )
            """
        )
        raw = cursor.execute(
            """
            SELECT
                authority,
                citation,
                species_taxon
            FROM gene
            GROUP BY
                authority,
                citation,
                species_taxon
            """
        ).fetchall()

        has_symbols = cursor.execute(
            """
            SELECT
                authority,
                citation,
                species_taxon,
                COUNT(symbol)
            FROM gene
            WHERE
                symbol IS NOT NULL
            GROUP BY
                authority,
                citation,
                species_taxon
            """
        ).fetchall()

        has_symbols = {
            row[:-1]: row[-1]
            for row in has_symbols
        }

        values = [
            (*row, 1)
            if row in has_symbols and has_symbols[row] > 0
            else (*row, 0)
            for row in raw
        ]

        cursor.executemany(
            f"""
            INSERT INTO {table_name} (
                authority,
                citation,
                species_taxon,
                has_symbols
            )
            VALUES (?, ?, ?, ?)
            """,
            values
        )
        db_utils.create_index(
            cursor=cursor,
            idx_name=index_name,
            table_name=table_name,
            column_tuple=("species_taxon", "authority", "has_symbols")
        )


def apply_mapping(
        gene_list,
        mapping,
        assign_placeholders=True,
        placeholder_prefix=None):
    """
    Apply a mapping to a gene list, only keeping 1:1 matches.

    Parameters
    ----------
    gene_list:
        list of str. The input genes
    mapping:
        a dict mapping the genes in gene_list
        to the lists of matches produced by the MMCGeneMapper
    assign_placeholders:
        a boolean.
        If True, genes that have zero or many matches will be given
        placeholder names. If false, their original names will
        be passed through. Genes with degenerate mappings will
        always be given placeholder names.
    placeholder_prefix:
        an optional string to be added to the placeholder
        names given to unmappable genes to indicate at what
        point in the mapping they became unmappable

    Returns
    -------
    A dict
        'failure_log': tally of how many genes failed and why
        'gene_list': the mapped gene list

    Notes
    -----
    Genes that fail to map in a 1:1 way will be assigned
    a name like UNMAPPABLE_{REASON}_{ii}

    Genes that are degenerate (i.e. that map to the same
    identifier) will also be marked as unmappable
    """
    failure_log = {
        'zero matches': 0,
        'many matches': 0
    }

    new_gene_list = []
    for gene in gene_list:
        this = mapping[gene]
        if len(this) == 1:
            assn = this[0]
        else:
            if len(this) == 0:
                ct = failure_log['zero matches']
                if assign_placeholders:
                    if placeholder_prefix is None:
                        assn = f'UNMAPPABLE_NO_MATCH_{ct}'
                    else:
                        assn = f'{placeholder_prefix}:UNMAPPABLE_NO_MATCH_{ct}'
                else:
                    assn = gene
                failure_log['zero matches'] += 1
            else:
                ct = failure_log['many matches']
                if assign_placeholders:
                    if placeholder_prefix is None:
                        assn = f'UNMAPPABLE_MANY_MATCHES_{ct}'
                    else:
                        assn = (
                            f'{placeholder_prefix}:UNMAPPABLE_MANY_MATCHES'
                            f'_{ct}'
                        )
                else:
                    assn = gene
                failure_log['many matches'] += 1
        new_gene_list.append(assn)

    (new_gene_list,
     n_degen) = mask_degenerate_genes(
                    new_gene_list,
                    placeholder_prefix=placeholder_prefix)

    failure_log['degenerate matches'] = n_degen

    return {
        'failure_log': failure_log,
        'gene_list': new_gene_list
    }


def mask_degenerate_genes(
        gene_list,
        placeholder_prefix=None):
    """
    Take a list of gene identifiers and replace any genes
    that have identical identifiers with unique placeholder
    identifiers.

    Parameters
    ----------
    gene_list:
        list of str. The input genes
    placeholder_prefix:
        an optional string to be added to the placeholder
        names given to degenerate genes

    Returns
    -------
    List of unique gene identifiers

    Number of degenerate genes found.
    """
    # make sure genes have unique names
    salt = 0
    unq, ct = np.unique(gene_list, return_counts=True)
    degen = (ct > 1)
    unq = unq[degen]
    degen_gene_to_idx = {
        u: ii for ii, u in enumerate(unq)
    }
    degen_gene_to_ct = {
        g: 0 for g in degen_gene_to_idx
    }
    n_degen = 0
    if len(degen_gene_to_idx) > 0:
        new_gene_list = []
        n_degen = len(degen_gene_to_idx)*2
        for ii in range(len(gene_list)):
            gene = gene_list[ii]
            if gene in degen_gene_to_idx:
                pair = degen_gene_to_idx[gene]
                salt = degen_gene_to_ct[gene]
                if placeholder_prefix is None:
                    assn = f'UNMAPPABLE_DEGENERATE_{pair}_{salt}'
                else:
                    assn = (
                        f'{placeholder_prefix}:'
                        f'UNMAPPABLE_DEGENERATE_{pair}_{salt}'
                    )
                new_gene_list.append(assn)
                degen_gene_to_ct[gene] += 1
            else:
                new_gene_list.append(gene)
    else:
        new_gene_list = copy.deepcopy(gene_list)
    return new_gene_list, n_degen


def detect_species_and_authority(
        db_path,
        gene_list):
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

    Returns
    --------
    A dict
        {
          "authority": np.array of strings indicating authority for genes
          "species": human readable name of species
          "species_taxon": integer indentifying species
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
            gene_list=gene_list[ensembl_idx]
        )
    else:
        ensembl_authority = None

    if len(ncbi_idx) > 0:
        ncbi_authority = _detect_species_and_authority(
            db_path=db_path,
            gene_list=gene_list[ncbi_idx]
        )
    else:
        ncbi_authority = None

    if ensembl_authority is None and ncbi_authority is None:
        species = None
        species_taxon = None
    elif ncbi_authority is None:
        species = ensembl_authority['species']
        species_taxon = ensembl_authority['species_taxon']
    elif ensembl_authority is None:
        species = ncbi_authority['species']
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
        species = ensembl_authority['species']
        species_taxon = ncbi_authority['species_taxon']

    return {
        "authority": authority,
        "species": species,
        "species_taxon": species_taxon
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
