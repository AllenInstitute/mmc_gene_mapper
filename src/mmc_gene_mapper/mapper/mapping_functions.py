"""
Functions that actually do the mapping of genes
"""
import json
import pathlib
import shutil
import sqlite3
import tempfile
import time

import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.utils.str_utils as str_utils
import mmc_gene_mapper.utils.typing_utils as typing_utils
import mmc_gene_mapper.metadata.classes as metadata_classes
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.download.download_manager as download_manager
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils
import mmc_gene_mapper.mapper.species_detection as species_utils
import mmc_gene_mapper.query_db.query as query_utils


def identifiers_from_symbols(
        db_path,
        gene_symbol_list,
        species,
        authority_name,
        assign_placeholders=True,
        placeholder_prefix=None):
    """
    Find the mapping that converts gene symbols into
    gene identifiers. Apply that mapping and return
    the list of relevant gene identifiers.

    Parameters
    ----------
    db_path:
        path to the database being queried
    gene_symbol_list:
        list of gene symbols
    species:
        Species representing the species we are working with
    authority_name:
        name of the authority in whose identifiers
        we want the genes listed
    assign_placeholders:
        a boolean. If True, assign placeholder names
        to any genes that cannot be mapped
    placeholder_prefix:
        optional prefix to apply to the placeholer names
        given to unmappable genes.

    Returns
    -------
    A dict
        {
          "metadata": {
              a dict describing the citation according
              to which these symbols map to these
              identifiers
          },
          "failure_log": {
             summary of how many genes failed to be mapped
             for what reasons
          }
          "gene_list": [
              list of mapped gene identifiers
          ]
        }
    """
    typing_utils.check_type(
        arg_name='identifiers_from_symbols.species',
        arg=species,
        expected_type=metadata_classes.Species
    )

    mapping = identifiers_from_symbols_mapping(
        db_path=db_path,
        gene_symbol_list=gene_symbol_list,
        species_name=species.name,
        authority_name=authority_name
    )

    mapped_result = mapper_utils.apply_mapping(
        gene_list=gene_symbol_list,
        mapping=mapping['mapping'],
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )

    result = {
        "metadata": mapping["metadata"],
        "failure_log": mapped_result["failure_log"],
        "gene_list": mapped_result["gene_list"]
    }
    return result



def identifiers_from_symbols_mapping(
        db_path,
        gene_symbol_list,
        species_name,
        authority_name):
    """
    Find the mapping that converts gene symbols into
    gene identifiers

    Parameters
    ----------
    db_path:
        path to the database being queries
    gene_symbol_list:
        list of gene symbols
    species_name:
        name of the species we are working with
    authority_name:
        name of the authority in whose identifiers
        we want the genes listed

    Returns
    -------
    A dict
        {
          "metadata": {
              a dict describing the citation according
              to which these symbols map to these
              identifiers
          },
          "mapping": {
              a dict keyed on input symbols. Each symbol
              maps to the list of all gene identifiers
              that are associated with that symbol according
              to the source described in "metadata"
          }
        }
    """
    query_utils.does_path_exist(db_path)

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        species_taxon = query_utils._get_species_taxon(
            cursor=cursor,
            species_name=species_name)

    result = query_utils.translate_gene_identifiers(
        db_path=db_path,
        src_column="symbol",
        dst_column="identifier",
        src_list=gene_symbol_list,
        authority_name=authority_name,
        species_taxon=species_taxon,
        chunk_size=500
    )
    return result


def equivalent_genes(
        db_path,
        input_authority,
        output_authority,
        gene_list,
        species_name,
        citation_name,
        assign_placeholders=True,
        placeholder_prefix=None):
    """
    Return a mapping between gene identifiers from
    different authorities (NCBI vs ENSEMBL)

    Parameters
    ----------
    db_path:
        the path to the database being queried
    input_authority:
        a str; the name of the authority in which the input
        genes are identified
    output_authority:
        a str; the name of the authority you want to convert
        the identifiers to
    gene_list:
        list of gene identifiers (in input_authority) to be
        mapped
    species_name:
        name of species we are working with
    citation_name:
        name of citation to use to assess gene eqivalence
    assign_placeholders:
        a boolean. If True, assign placeholder names
        to any genes that cannot be mapped
    placeholder_prefix:
        optional prefix to apply to the placeholer names
        given to unmappable genes.

    Returns
    -------
    A dict
        {
          "metadata": {
              a dict describing the citation according
              to which these symbols map to these
              identifiers
          },
          "failure_log": {
             summary of how many genes failed to be mapped
             for what reasons
          }
          "gene_list": [
              list of mapped gene identifiers
          ]
        }
    """

    mapping = equivalent_genes_mapping(
        db_path=db_path,
        input_authority=input_authority,
        output_authority=output_authority,
        gene_list=gene_list,
        species_name=species_name,
        citation_name=citation_name
    )

    mapped_result = mapper_utils.apply_mapping(
        gene_list=gene_list,
        mapping=mapping['mapping'],
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )

    result = {
        'metadata': mapping['metadata'],
        'failure_log': mapped_result['failure_log'],
        'gene_list': mapped_result['gene_list']
    }

    return result


def ortholog_genes_mapping(
        db_path,
        authority,
        src_species,
        dst_species,
        gene_list,
        citation_name):
    """
    Return a mapping between gene identifiers from
    different species

    Parameters
    ----------
    db_path:
        path to the database being queried
    authority:
        a str; the name of the authority (ENSEMBL, NCBI etc.)
        we are working in
    src_species:
        a Species; the species we are starting from
    dst_species:
        as Species; the species we are mapping to
    gene_list:
        list of gene identifiers (in src_species) to be
        mapped
    citation_name:
        name of citation to use to assess gene eqivalence

    Returns
    -------
    A dict
        {
          "metadata": {
              a dict describing the citation according
              to which these symbols map to these
              identifiers
          },
          "mapping": {
              a dict keyed on input genes. Each gene
              maps to the list of all genes
              that are considered orthologs to that
              gene according to the source described
              in "metadata"
          }
        }
    """
    typing_utils.check_many_types(
        name_arr=['ortholog_genes_mapping.src_species',
                  'ortholog_genes_mapping.dst_species'],
        arg_arr=[src_species, dst_species],
        type_arr=[metadata_classes.Species, metadata_classes.Species]
    )

    return query_utils.get_ortholog_genes_from_identifiers(
        db_path=db_path,
        authority_name=authority,
        src_species_taxon=src_species.taxon,
        dst_species_taxon=dst_species.taxon,
        src_gene_list=gene_list,
        citation_name=citation_name,
        chunk_size=500
    )


def ortholog_genes(
        db_path,
        authority,
        src_species,
        dst_species,
        gene_list,
        citation_name,
        assign_placeholders=True,
        placeholder_prefix=None):
    """
    Return a mapping between gene identifiers from
    different species

    Parameters
    ----------
    db_path:
        the path to the database being queried
    authority:
        a str; the name of the authority (ENSEMBL, NCBI etc.)
        we are working in
    src_species:
        a Species; the name of the species we are starting from
    dst_species:
        as Species; the name of the species we are mapping to
    gene_list:
        list of gene identifiers (in src_species) to be
        mapped
    citation_name:
        name of citation to use to assess gene eqivalence
    assign_placeholders:
        a boolean. If True, assign placeholder names
        to any genes that cannot be mapped
    placeholder_prefix:
        optional prefix to apply to the placeholer names
        given to unmappable genes.

    Returns
    -------
    A dict
        {
          "metadata": {
              a dict describing the citation according
              to which these symbols map to these
              identifiers
          },
          "failure_log": {
             summary of how many genes failed to be mapped
             for what reasons
          }
          "gene_list": [
              list of mapped gene identifiers
          ]
        }
    """
    mapping = ortholog_genes_mapping(
        db_path=db_path,
        authority=authority,
        src_species=src_species,
        dst_species=dst_species,
        gene_list=gene_list,
        citation_name=citation_name
    )

    mapped_result = mapper_utils.apply_mapping(
        gene_list=gene_list,
        mapping=mapping['mapping'],
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )

    return {
        'metadata': mapping['metadata'],
        'failure_log': mapped_result['failure_log'],
        'gene_list': mapped_result['gene_list']
    }



def equivalent_genes_mapping(
        db_path,
        input_authority,
        output_authority,
        gene_list,
        species_name,
        citation_name):
    """
    Return a mapping between gene identifiers from
    different authorities (NCBI vs ENSEMBL)

    Parameters
    ----------
    db_path:
        path to the database being queried
    input_authority:
        a str; the name of the authority in which the input
        genes are identified
    output_authority:
        a str; the name of the authority you want to convert
        the identifiers to
    gene_list:
        list of gene identifiers (in input_authority) to be
        mapped
    species_name:
        name of species we are working with
    citation_name:
        name of citation to use to assess gene eqivalence

    Returns
    -------
    A dict
        {
          "metadata": {
              a dict describing the citation according
              to which these symbols map to these
              identifiers
          },
          "mapping": {
              a dict keyed on input genes. Each gene
              maps to the list of all genes
              that are considered equivalent to that
              gene according to the source described
              in "metadata"
          }
        }
    """
    return query_utils.get_equivalent_genes_from_identifiers(
        db_path=db_path,
        input_authority_name=input_authority,
        output_authority_name=output_authority,
        input_gene_list=gene_list,
        species_name=species_name,
        citation_name=citation_name,
        chunk_size=500
    )
