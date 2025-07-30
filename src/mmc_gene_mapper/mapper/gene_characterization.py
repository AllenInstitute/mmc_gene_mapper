import numpy as np
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.utils.str_utils as str_utils


def characterize_gene_identifiers(
        db_path,
        gene_list,
        chunk_size=10000):
    """
    Take a list of genes identifiers. Determine which are symbols,
    which are ENSEMBL IDs and which are NCBI IDs.

    Parameters
    ----------
    db_path:
        path to database backing gene mapper
    gene_list:
        list of gene identifiers being characterized
    chunk_size:
        number of genes to check at a time

    Returns
    -------
    numpy array of strings characterizing genes as
    'ENSEMBL', 'NCBI' or 'symbol'

    Notes
    -----
    The function will first try to characterize the genes by comparing
    them to known gene identifiers and gene symbols. Anything remaining
    will be characterized according to simple regex matching
    """
    file_utils.assert_is_file(db_path)
    gene_list = np.array(gene_list)
    result = np.array([None]*len(gene_list))
    n_genes = len(gene_list)
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        for i0 in range(0, n_genes, chunk_size):
            i1 = min(n_genes, i0+chunk_size)
            chunk = gene_list[i0: i1]

            # see which genes in gene_list are known identifiers
            identifier_to_authority = _get_authority_lookup(
                cursor=cursor,
                chunk=chunk
            )

            chunk_idx = np.array(
                [ii-i0 for ii in range(i0, i1, 1)
                 if result[ii] is None]
            )
            remaining_chunk = chunk[chunk_idx]

            # see which of remaining are known symbols
            symbol_list = _get_symbol_list(
                cursor=cursor,
                chunk=remaining_chunk)

            for ii, gene in enumerate(chunk):
                if gene in identifier_to_authority:
                    result[ii+i0] = identifier_to_authority[gene]
                elif gene in symbol_list:
                    result[ii+i0] = "symbol"

    unknown_idx = np.array(
        [ii for ii in range(n_genes) if result[ii] is None]
    )
    if len(unknown_idx) > 0:
        result[unknown_idx] = str_utils.characterize_gene_identifiers_by_re(
            gene_list[unknown_idx]
        )
    return result


def _get_authority_lookup(
        cursor,
        chunk):
    """
    Query gene table on identifier column for values in
    chunk. Return a lookup table mapping values in chunk
    to authority. If there are degeneracies, raise an
    exception.
    """
    query = """
        SELECT
            gene.identifier,
            authority.name
        FROM
            gene
        JOIN authority ON
            gene.authority = authority.id
        WHERE
            gene.identifier IN (
    """
    query += ", ".join(["?"]*len(chunk))
    query += ")"
    result = cursor.execute(
        query,
        chunk
    ).fetchall()
    lookup = dict()
    err_msg = ""
    for row in result:
        gene = row[0]
        authority = row[1]
        if gene in lookup:
            if lookup[gene] != authority:
                err_msg += f"{gene}: {lookup[gene]} and {authority}\n"
        else:
            lookup[gene] = authority
    if len(err_msg) > 0:
        msg = "The following genes mapped to more than one authority\n"
        msg += f"{err_msg}"
        raise MultipleAuthorityError(
            msg
        )
    return lookup


def _get_symbol_list(
        cursor,
        chunk):
    """
    Check which genes in chunk are known symbols in the gene table.
    Return list of valid symbos.
    """
    query = """
        SELECT
            DISTINCT(symbol)
        FROM
            gene
        WHERE
            symbol IN (
    """
    query += ", ".join(["?"]*len(chunk))
    query += ")"
    result = cursor.execute(query, chunk).fetchall()
    return [r[0] for r in result]


class MultipleAuthorityError(Exception):
    pass
