"""
Define functions to convert a set of genes of disparate
authority into a single authority (ENSEMBL or NCBI)


first step is to move mapping finders, etc. out of
the mapper class into a function module
"""

import numpy as np
import sqlite3

import mmc_gene_mapper.metadata.classes as metadata_classes
import mmc_gene_mapper.query_db.query as query_utils
import mmc_gene_mapper.mapper.species_detection as species_detection
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils
import mmc_gene_mapper.mapper.mapping_functions as mapping_functions


def arbitrary_mapping(
        db_path,
        gene_list,
        dst_species,
        dst_authority,
        ortholog_citation='NCBI'):
    """
    Parameters
    ----------
    db_path:
        path to the database being queries
    gene_list:
        list of gene identifiers being mapped
    dst_species:
        name of species being mapped to
    dst_authority:
        name of authority being mapped to
    ortholog_citation:
        citation to use for ortholog mapping, if necessary
    """
    query_utils.does_path_exist(db_path)

    if dst_authority not in ('NCBI', 'ENSEMBL'):
        raise ValueError(
            f"Unclear how to map to authority '{dst_authority}'; "
            "Must be either 'NCBI' or 'ENSEMBL'"
        )

    n_genes = len(gene_list)
    metadata = []

    src_authority = species_detection.detect_species_and_authority(
        db_path=db_path,
        gene_list=gene_list
    )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        dst_species_taxon = query_utils._get_species_taxon(
            cursor=cursor,
            species_name=dst_species
        )

    if src_authority['species'] is None:
        src_authority['species'] = dst_species
        src_authority['species_taxon'] = dst_species_taxon

    src_species_taxon = src_authority['species_taxon']

    if src_species_taxon is None or dst_species_taxon == src_species_taxon:
        need_orthologs = False
    else:
        need_orthologs = True

    if need_orthologs:
        current = convert_authority_in_bulk(
            db_path=db_path,
            gene_list=gene_list,
            src_authority=src_authority,
            dst_authority='NCBI'
        )

        for auth_metadata in current['metadata']:
            metadata.append(auth_metadata)

        current = mapping_functions.ortholog_genes(
            db_path=db_path,
            authority='NCBI',
            src_species_name=src_species_taxon,
            dst_species_name=dst_species_taxon,
            gene_list=current['gene_list'],
            citation_name=ortholog_citation,
            assign_placeholders=True,
            placeholder_prefix='ortholog'
        )

        current_authority = {
            'authority': np.array(['NCBI']*n_genes),
            'species': dst_species,
            'species_taxon': dst_species_taxon
        }

        current['metadata']['mapping'] = {
            'axis': 'species',
            'from': {
                'name': src_authority['species'],
                'taxon': src_authority['species_taxon']
            },
            'to': {
                'name': dst_species,
                'taxon': dst_species_taxon
            }
        }

        metadata.append(current['metadata'])

    else:
        current_authority = src_authority
        current = {
            'gene_list': gene_list,
        }

    if set(current_authority['authority']) != set([dst_authority]):
        current = convert_authority_in_bulk(
            db_path=db_path,
            gene_list=current['gene_list'],
            src_authority=current_authority,
            dst_authority=dst_authority
        )

        for auth_metadata in current['metadata']:
            metadata.append(auth_metadata)

    return {
        'metadata': metadata,
        'gene_list': current['gene_list']
    }


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

    metadata = []
    failure_log = None

    if len(symbol_idx) > 0:
        raw = mapping_functions.identifiers_from_symbols(
            db_path=db_path,
            gene_symbol_list=gene_list[symbol_idx],
            species_name=species_taxon,
            authority_name=dst_authority,
            assign_placeholders=True,
            placeholder_prefix=f"symbol:{dst_authority}"
        )
        result[symbol_idx] = np.array(raw['gene_list'])
        metadata.append(raw['metadata'])
        failure_log = raw['failure_log']

    for input_authority, idx_arr in [('ENSEMBL', ensembl_idx),
                                     ('NCBI', ncbi_idx)]:
        if len(idx_arr) == 0:
            continue

        if input_authority == dst_authority:
            result[idx_arr] = gene_list[idx_arr]
        else:
            raw = mapping_functions.equivalent_genes(
                db_path=db_path,
                input_authority=input_authority,
                output_authority=dst_authority,
                gene_list=gene_list[idx_arr],
                species_name=species_taxon,
                citation_name='NCBI',
                assign_placeholders=True,
                placeholder_prefix=f"{input_authority}:{dst_authority}"
            )
            result[idx_arr] = np.array(raw['gene_list'])
            metadata.append(raw['metadata'])
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

    if failure_log is not None:
        failure_log['degenerate matches'] += n_degen
    else:
        failure_log = {'degenerate matches': n_degen}

    return {
        'metadata': metadata,
        'failure_log': failure_log,
        'gene_list': list(result)
    }
