"""
Define functions to convert a set of genes of disparate
authority into a single authority (ENSEMBL or NCBI)


first step is to move mapping finders, etc. out of
the mapper class into a function module
"""

import numpy as np

import mmc_gene_mapper.mapper.mapper_utils as mapper_utils
import mmc_gene_mapper.mapper.mapping_functions as mapping_functions


def convert_authority_in_bulk(
        db_path,
        gene_list,
        src_authority,
        dst_authority):
    """
    db_path:
        path to database being queried
    gene_list:
        genes being mapped
    src_authority:
        dict containing species and authority
        information for the genes (output of
        detect_species_and_authority)
    dst_authority:
        authority to which genes are being mapped
    """
    if dst_authority not in ('NCBI', 'ENSEMBL'):
        raise ValueError(
            f"Unclear how to map to authority '{dst_authority}'; "
            "Must be either 'NCBI' or 'ENSEMBL'"
        )
    gene_list = np.array(gene_list)

    symbol_idx = np.where(
        src_authority['authority'] == 'symbol'
    )[0]
    ensembl_idx = np.where(
        src_authority['authority'] == 'ENSEMBL'
    )[0]
    ncbi_idx = np.where(
        src_authority['authority'] == 'NCBI'
    )[0]
    n_genes = len(gene_list)

    n_found = (
        len(symbol_idx)
        + len(ensembl_idx)
        + len(ncbi_idx)
    )

    if n_found != n_genes:
        raise RuntimeError(
            "Not all genes accounted for in authorities "
            "('symbol', 'ENSEMBL', 'NCBI'); "
            f"input n: {n_genes}\n"
            f"n authorities: {n_found}"
        )

    result = np.array([None]*n_genes)
    species_taxon = src_authority['species_taxon']

    metadata = dict()
    failure_log = None

    if len(symbol_idx) > 0:
        prefix=f'symbol:{dst_authority}'
        raw = mapping_functions.identifiers_from_symbols(
            db_path=db_path,
            gene_symbol_list=gene_list[symbol_idx],
            species_name=species_taxon,
            authority_name=dst_authority,
            assign_placeholders=True,
            placeholder_prefix=prefix
        )
        result[symbol_idx] = np.array(raw['gene_list'])
        metadata[prefix] = raw['metadata']
        failure_log = raw['failure_log']

    for input_authority, idx_arr in [('ENSEMBL', ensembl_idx),
                                     ('NCBI', ncbi_idx)]:
        if len(idx_arr) == 0:
            continue

        if input_authority == dst_authority:
            result[idx_arr] = gene_list[idx_arr]
        else:
            prefix = f'{input_authority}:{dst_authority}'
            raw = mapping_functions.equivalent_genes(
                db_path=db_path,
                input_authority=input_authority,
                output_authority=dst_authority,
                gene_list=gene_list[idx_arr],
                species_name=species_taxon,
                citation_name='NCBI',
                assign_placeholders=True,
                placeholder_prefix=prefix
            )
            result[idx_arr] = np.array(raw['gene_list'])
            metadata[prefix] = raw['metadata']
            if failure_log is None:
                failure_log = raw['failure_log']
            else:
                for key in failure_log:
                    failure_log[key] += raw['failure_log'][key]

    (result,
     n_degen) = mapper_utils.mask_degenerate_genes(
        gene_list=result,
        placeholder_prefix=f'{dst_authority}'
    )

    failure_log['degenerate matches'] += n_degen

    return {
        'metadata': metadata,
        'failure_log': failure_log,
        'gene_list': list(result)
    }
    

            
