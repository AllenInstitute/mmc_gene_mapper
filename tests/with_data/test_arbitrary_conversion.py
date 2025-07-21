"""
Test our "as generic as possible" gene mapping functions
"""
import pytest

import itertools
import numpy as np
import sqlite3

import mmc_gene_mapper.metadata.classes as metadata_classes
import mmc_gene_mapper.query_db.query as query_utils
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

    src_gene_data = species_detection.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=gene_list
    )

    result = arbitrary_conversion._convert_authority_in_bulk(
        db_path=mapper_fixture.db_path,
        gene_list=gene_list,
        src_gene_data=src_gene_data,
        dst_authority=dst_authority
    )

    np.testing.assert_array_equal(
        np.array(expected),
        np.array(result['gene_list'])
    )


@pytest.mark.parametrize(
    "gene_list, dst_authority, expected, use_class",
    [
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22"],
       "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11"],
       False),
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22",
        "symbol:777", "ENSX10009", "NCBIGene:671321"],
       "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11",
        "symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
        "ENSEMBL:NCBI:UNMAPPABLE_NO_MATCH_0",
        "NCBIGene:671321"],
       False),
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "symbol:13",
        "ENSX22"],
       "NCBI",
       ["NCBIGene:12", "NCBIGene:14",
        "NCBI:UNMAPPABLE_DEGENERATE_0_0",
        "NCBIGene:16",
        "NCBI:UNMAPPABLE_DEGENERATE_0_1",
        "NCBIGene:11"],
       False),
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22"],
       "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11"],
       True),
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "ENSX22",
        "symbol:777", "ENSX10009", "NCBIGene:671321"],
       "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13", "NCBIGene:16",
        "NCBIGene:11",
        "symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
        "ENSEMBL:NCBI:UNMAPPABLE_NO_MATCH_0",
        "NCBIGene:671321"],
       True),
      (["symbol:12", "symbol:14", "ENSX26", "NCBIGene:16", "symbol:13",
        "ENSX22"],
       "NCBI",
       ["NCBIGene:12", "NCBIGene:14",
        "NCBI:UNMAPPABLE_DEGENERATE_0_0",
        "NCBIGene:16",
        "NCBI:UNMAPPABLE_DEGENERATE_0_1",
        "NCBIGene:11"],
       True
       ),
      (["symbol:12", "symbol:14", "symbol:13"],
       "NCBI",
       ["NCBIGene:12", "NCBIGene:14", "NCBIGene:13"],
       True
       ),
      (["symbol:22", "symbol:26", "symbol:24", "symbol:28"],
       "ENSEMBL",
       ["symbol:ENSEMBL:UNMAPPABLE_NO_MATCH_0",
        "symbol:ENSEMBL:UNMAPPABLE_NO_MATCH_1",
        "ENSX26",
        "ENSX30"],
       True
       ),
    ]
)
def test_arbitrary_mapping_function_no_ortholog(
        gene_list,
        dst_authority,
        expected,
        mapper_fixture,
        use_class):

    if use_class:
        result = mapper_fixture.map_genes(
            gene_list=gene_list,
            dst_species='jabberwock',
            dst_authority=dst_authority,
            ortholog_citation='NCBI'
        )
    else:
        with sqlite3.connect(mapper_fixture.db_path) as conn:
            cursor = conn.cursor()
            dst_species = query_utils.get_species(
                cursor=cursor,
                species='jabberwock'
            )

        assert isinstance(dst_species, metadata_classes.Species)

        result = arbitrary_conversion.arbitrary_mapping(
            db_path=mapper_fixture.db_path,
            gene_list=gene_list,
            dst_species=dst_species,
            dst_authority=metadata_classes.Authority(dst_authority),
            ortholog_citation='NCBI'
        )

    np.testing.assert_array_equal(
        np.array(result['gene_list']),
        np.array(expected)
    )


@pytest.mark.parametrize(
    "gene_list, dst_authority, expected, use_class",
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
        "ortholog:UNMAPPABLE_NO_MATCH_1"],
       False),
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
        "ENSX14"],
       False),
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
        "symbol:NCBI:UNMAPPABLE_NO_MATCH_0"],
       False),
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
        "ortholog:UNMAPPABLE_NO_MATCH_1"],
       True),
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
        "ENSX14"],
       True),
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
        "symbol:NCBI:UNMAPPABLE_NO_MATCH_0"],
       True),
    ]
)
def test_arbitrary_mapping_function_yes_ortholog(
        gene_list,
        dst_authority,
        expected,
        mapper_fixture,
        use_class):

    if use_class:
        result = mapper_fixture.map_genes(
            gene_list=gene_list,
            dst_species='human',
            dst_authority=dst_authority,
            ortholog_citation='NCBI'
        )
    else:
        with sqlite3.connect(mapper_fixture.db_path) as conn:
            cursor = conn.cursor()
            dst_species = query_utils.get_species(
                cursor=cursor,
                species='human'
            )
        result = arbitrary_conversion.arbitrary_mapping(
            db_path=mapper_fixture.db_path,
            gene_list=gene_list,
            dst_species=dst_species,
            dst_authority=metadata_classes.Authority(dst_authority),
            ortholog_citation='NCBI'
        )

    np.testing.assert_array_equal(
        np.array(result['gene_list']),
        np.array(expected)
    )


@pytest.mark.parametrize(
    "gene_list, dst_authority, use_class",
    itertools.product(
        (["symbol:12", "symbol:14", "ENSX26",
          "NCBIGene:16", "ENSX22"],
         ["symbol:12", "symbol:14", "symbol:13"]),
        ("NCBI", "ENSEMBL"),
        (True, False)
    )
)
def test_arbitrary_mapping_function_with_log(
        gene_list,
        dst_authority,
        mapper_fixture,
        use_class):
    """
    Make sure that mapper can accept a logger class.
    This is just a smoke test.
    """

    if "NCBIGene:16" in gene_list:
        dst_species_name = "human"
    else:
        dst_species_name = "jabberwock"

    class DummyLog(object):
        def info(self, msg, to_stdout=False):
            if to_stdout:
                print(msg)
            pass

    if use_class:
        mapper_fixture.map_genes(
            gene_list=gene_list,
            dst_species=dst_species_name,
            dst_authority=dst_authority,
            ortholog_citation='NCBI',
            log=DummyLog()
        )
    else:
        with sqlite3.connect(mapper_fixture.db_path) as conn:
            cursor = conn.cursor()
            dst_species = query_utils.get_species(
                cursor=cursor,
                species=dst_species_name
            )
        arbitrary_conversion.arbitrary_mapping(
            db_path=mapper_fixture.db_path,
            gene_list=gene_list,
            dst_species=dst_species,
            dst_authority=metadata_classes.Authority(dst_authority),
            ortholog_citation='NCBI',
            log=DummyLog()
        )


def test_arbitrary_conversion_typing():

    msg = "dst_species must be of type"
    with pytest.raises(ValueError, match=msg):
        arbitrary_conversion.arbitrary_mapping(
            db_path='not/a/file.db',
            gene_list=[],
            dst_species='human',
            dst_authority=metadata_classes.Authority('NCBI')
        )

    msg = "dst_authority must be of type"
    with pytest.raises(ValueError, match=msg):
        arbitrary_conversion.arbitrary_mapping(
            db_path='not/a/file.db',
            gene_list=[],
            dst_species=metadata_classes.Species(name='blerg', taxon=1),
            dst_authority='NCBI'
        )
