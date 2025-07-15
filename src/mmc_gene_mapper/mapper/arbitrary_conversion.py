"""
Define functions to convert a set of genes of disparate
authority into a single authority (ENSEMBL or NCBI)


first step is to move mapping finders, etc. out of
the mapper class into a function module
"""

import numpy as np
import sqlite3

import mmc_gene_mapper.utils.typing_utils as typing_utils
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
        metadata_classes.Species representing
        the species to which we are mapping
    dst_authority:
        metadata_classes.Authority representing
        the desired output authority
    ortholog_citation:
        citation to use for ortholog mapping, if necessary
    """
    typing_utils.check_many_types(
        name_arr=["dst_species", "dst_authority"],
        arg_arr=[dst_species, dst_authority],
        type_arr=[metadata_classes.Species,
                  metadata_classes.Authority]
    )
    query_utils.does_path_exist(db_path)

    if dst_authority.name not in ('NCBI', 'ENSEMBL'):
        raise ValueError(
            f"Unclear how to map to authority '{dst_authority}'; "
            "Must be either 'NCBI' or 'ENSEMBL'"
        )

    n_genes = len(gene_list)
    metadata = []

    src_gene_data = species_detection.detect_species_and_authority(
        db_path=db_path,
        gene_list=gene_list
    )

    print('====got src_gene_data====')

    if src_gene_data['species'] is None:
        src_gene_data['species'] = dst_species

    src_species = src_gene_data['species']

    if dst_species.taxon == src_species.taxon:
        need_orthologs = False
    else:
        need_orthologs = True

    print(f'====need_orthologs {need_orthologs}====')

    if need_orthologs:
        current = _convert_authority_in_bulk(
            db_path=db_path,
            gene_list=gene_list,
            src_gene_data=src_gene_data,
            dst_authority='NCBI'
        )

        print('=====converted authority before orthologs=====')

        for auth_metadata in current['metadata']:
            metadata.append(auth_metadata)

        current = mapping_functions.ortholog_genes(
            db_path=db_path,
            authority='NCBI',
            src_species=src_species,
            dst_species=dst_species,
            gene_list=current['gene_list'],
            citation_name=ortholog_citation,
            assign_placeholders=True,
            placeholder_prefix='ortholog'
        )

        current_gene_data = {
            'authority': np.array(['NCBI']*n_genes),
            'species': dst_species
        }

        metadata.append(current['metadata'])
        print('======mapped orthologs======')

    else:
        current_gene_data = src_gene_data
        current = {
            'gene_list': gene_list,
        }

    if set(current_gene_data['authority']) != set([dst_authority]):
        print('=======performing final authority conversion=======')
        current = _convert_authority_in_bulk(
            db_path=db_path,
            gene_list=current['gene_list'],
            src_gene_data=current_gene_data,
            dst_authority=dst_authority.name
        )
        print('=======done with final authority conversion======')

        for auth_metadata in current['metadata']:
            metadata.append(auth_metadata)

    return {
        'metadata': metadata,
        'gene_list': current['gene_list']
    }


def _convert_authority_in_bulk(
        db_path,
        gene_list,
        src_gene_data,
        dst_authority):
    """
    db_path:
        path to database being queried
    gene_list:
        genes being mapped
    src_gene_data:
        dict containing species and authority
        information for the genes. Specifically
            {'authority': [list of strings],
             'species': an instance of Species
            }
    dst_authority:
        authority to which genes are being mapped
    """
    species = src_gene_data['species']

    if dst_authority not in ('NCBI', 'ENSEMBL'):
        raise ValueError(
            f"Unclear how to map to authority '{dst_authority}'; "
            "Must be either 'NCBI' or 'ENSEMBL'"
        )
    gene_list = np.array(gene_list)

    symbol_idx = np.where(
        src_gene_data['authority'] == 'symbol'
    )[0]
    ensembl_idx = np.where(
        src_gene_data['authority'] == 'ENSEMBL'
    )[0]
    ncbi_idx = np.where(
        src_gene_data['authority'] == 'NCBI'
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

    metadata = []
    failure_log = None

    if len(symbol_idx) > 0:
        raw = mapping_functions.identifiers_from_symbols(
            db_path=db_path,
            gene_symbol_list=gene_list[symbol_idx],
            species=species,
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
                species=species,
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
