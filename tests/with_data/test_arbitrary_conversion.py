"""
Test our "as generic as possible" gene mapping functions
"""
import pytest

import numpy as np

import mmc_gene_mapper.mapper.species_detection as species_detection
import mmc_gene_mapper.mapper.arbitrary_conversion as arbitrary_conversion


@pytest.mark.parametrize(
    "gene_list, dst_authority, expected",
    [
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22"],
        "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11"]),
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22",
        "symbol:777", "ENSX10009", "NCBIGene:671321"],
        "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11",
        "symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
        "ENSEMBL:NCBI:UNMAPPABLE_NO_MATCH_0",
        "NCBIGene:671321"]),
       (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "symbol:13",
         "ENSX22"],
        "NCBI",
       ["NCBIGene:12", "NCBIGene:14",
        "NCBI:UNMAPPABLE_DEGENERATE_0_0",
        "NCBIGene:16",
        "NCBI:UNMAPPABLE_DEGENERATE_0_1",
        "NCBIGene:11"]),
    ]
)
def test_convert_authority_in_bulk(
        gene_list,
        dst_authority,
        expected,
        mapper_fixture):

    src_authority = species_detection.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=gene_list
    )

    result = arbitrary_conversion.convert_authority_in_bulk(
        db_path=mapper_fixture.db_path,
        gene_list=gene_list,
        src_authority=src_authority,
        dst_authority=dst_authority
    )

    np.testing.assert_array_equal(
        np.array(expected),
        np.array(result['gene_list'])
    )


@pytest.mark.parametrize(
    "gene_list, dst_authority, expected",
    [
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22"],
        "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11"]),
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22",
        "symbol:777", "ENSX10009", "NCBIGene:671321"],
        "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11",
        "symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
        "ENSEMBL:NCBI:UNMAPPABLE_NO_MATCH_0",
        "NCBIGene:671321"]),
       (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "symbol:13",
         "ENSX22"],
        "NCBI",
       ["NCBIGene:12", "NCBIGene:14",
        "NCBI:UNMAPPABLE_DEGENERATE_0_0",
        "NCBIGene:16",
        "NCBI:UNMAPPABLE_DEGENERATE_0_1",
        "NCBIGene:11"]),
    ]
)
def test_arbitrary_mapping_function_no_ortholog(
        gene_list,
        dst_authority,
        expected,
        mapper_fixture):

    result = arbitrary_conversion.arbitrary_mapping(
        db_path=mapper_fixture.db_path,
        gene_list=gene_list,
        dst_species='jabberwock',
        dst_authority=dst_authority,
        ortholog_citation='NCBI'
    )

    np.testing.assert_array_equal(
        np.array(result['gene_list']),
        np.array(expected)
    )


@pytest.mark.parametrize(
    "gene_list, dst_authority, expected",
    [
      (["symbol:12",
        "symbol:14",
        "ENSX26",
        "NCBIGene:16",
        "ENSX22"],
        "NCBI",
       ["NCBIGene:4",
        "NCBIGene:0",
        "NCBIGene:1",
        "ortholog:UNMAPPABLE_NO_MATCH_0",
        "ortholog:UNMAPPABLE_NO_MATCH_1"]),
      (["symbol:12",
        "symbol:14",
        "ENSX26",
        "NCBIGene:16",
        "ENSX22",
        "NCBIGene:15"],
        "ENSEMBL",
       ["NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_0",
        "NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_1",
        "ENSX2",
        "ortholog:UNMAPPABLE_NO_MATCH_0",
        "ortholog:UNMAPPABLE_NO_MATCH_1",
        "ENSX14"]),
      (["symbol:12",
        "symbol:14",
        "ENSX26",
        "NCBIGene:16",
        "ENSX22",
        "symbol:15",
        "ENSX88888",
        "symbol:987654"],
        "ENSEMBL",
       ["NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_0",
        "NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_1",
        "ENSX2",
        "ortholog:UNMAPPABLE_NO_MATCH_0",
        "ortholog:UNMAPPABLE_NO_MATCH_1",
        "ENSX14",
        "ENSEMBL:NCBI:UNMAPPABLE_NO_MATCH_0",
        "symbol:NCBI:UNMAPPABLE_NO_MATCH_0"]),
    ]
)
def test_arbitrary_mapping_function_yes_ortholog(
        gene_list,
        dst_authority,
        expected,
        mapper_fixture):

    result = arbitrary_conversion.arbitrary_mapping(
        db_path=mapper_fixture.db_path,
        gene_list=gene_list,
        dst_species='human',
        dst_authority=dst_authority,
        ortholog_citation='NCBI'
    )

    np.testing.assert_array_equal(
        np.array(result['gene_list']),
        np.array(expected)
    )
