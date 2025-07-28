import pytest

import numpy as np

import mmc_gene_mapper.mapper.species_detection as species_detection
import mmc_gene_mapper.mapper.mapper as mapper_module


@pytest.mark.parametrize(
    "gene_list, dst_species, dst_authority, expected, err_flavor, err_msg",
    [
     # within human NCBI->ENSEMBL
     (["NCBIGene:100", "NCBIGene:112", "NCBIGene:113", "NCBIGene:124"],
      "human",
      "ENSEMBL",
      ["ENSG100", "ENSG112", "NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_0", "ENSG124"],
      None,
      None),
     # within human NCBI->ENSEMBL, with duplicate gene
     (["NCBIGene:136", "NCBIGene:100", "NCBIGene:112", "NCBIGene:113",
       "NCBIGene:124", "NCBIGene:100", "NCBIGene:136"],
      "human",
      "ENSEMBL",
      ["NCBI:ENSEMBL:UNMAPPABLE_DEGENERATE_1_0",
       "NCBI:ENSEMBL:UNMAPPABLE_DEGENERATE_0_0",
       "ENSG112",
       "NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_0",
       "ENSG124",
       "NCBI:ENSEMBL:UNMAPPABLE_DEGENERATE_0_1",
       "NCBI:ENSEMBL:UNMAPPABLE_DEGENERATE_1_1"],
      None,
      None),
     # within human, symbol -> ENSEMBL; will error because src could be mouse
     (["symbol_12", "symbol_24", "symbol_74"],
      "human",
      "ENSEMBL",
      None,
      species_detection.InconsistentSpeciesError,
      "The gene symbols you gave are consistent with "
      "more than one species"),
     # human -> mouse based just on symbols
     # recall can only ortholog through fish if i%2 == i%3 == 0
     (["NCBIGene:999999", "symbol_12", "symbol_24",
       "symbol_74", "symbol_8888"],
      "mouse",
      "ENSEMBL",
      ["symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
       "ENSM12",
       "ENSM24",
       "ortholog:UNMAPPABLE_NO_MATCH_1",
       "ortholog:UNMAPPABLE_NO_MATCH_2"],
      None,
      None),
     # human -> mouse based just on symbols
     # recall can only ortholog through fish if i%2 == i%3 == 0
     (["NCBIGene:999999",
       "symbol_12",
       "symbol_24",
       "symbol_74",
       "symbol_6",
       "symbol_8888"],
      "mouse",
      "NCBI",
      ["symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
       "NCBIGene:12",
       "NCBIGene:24",
       "ortholog:UNMAPPABLE_NO_MATCH_1",
       "NCBIGene:6",
       "ortholog:UNMAPPABLE_NO_MATCH_2"],
      None,
      None),
     # human -> mouse based just on symbols;
     # use ENSEMBL-specific symbol, which will not map to NCBI
     (["NCBIGene:999999",
       "symbol_12",
       "e_human_24",
       "symbol_74",
       "symbol_6",
       "symbol_8888"],
      "mouse",
      "NCBI",
      ["symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
       "NCBIGene:12",
       "symbol:NCBI:UNMAPPABLE_NO_MATCH_1",
       "ortholog:UNMAPPABLE_NO_MATCH_2",
       "NCBIGene:6",
       "ortholog:UNMAPPABLE_NO_MATCH_3"],
      None,
      None),
     # within human, symbol->ENSEMBL using ensembl-only symbol;
     # (in case we ever allow daisy chain symbl->ENS->NCBI->orth...)
     (["NCBIGene:999999",
       "symbol_12",
       "e_human_24",
       "symbol_74",
       "symbol_6",
       "symbol_8888"],
      "human",
      "ENSEMBL",
      ["symbol:ENSEMBL:UNMAPPABLE_NO_MATCH_0",
       "ENSG112",
       "ENSG124",
       "ENSG174",
       "ENSG106",
       "ENSG8888"],
      None,
      None),
    ]
)
def test_arbitrary_mapping(
        mapper_db_path_fixture,
        gene_list,
        dst_species,
        dst_authority,
        expected,
        err_flavor,
        err_msg):

    mapper = mapper_module.MMCGeneMapper(
        db_path=mapper_db_path_fixture
    )

    if err_flavor is None:
        actual = mapper.map_genes(
            gene_list=gene_list,
            dst_species=dst_species,
            dst_authority=dst_authority,
            ortholog_citation='NCBI',
            log=None,
            invalid_mapping_prefix=None
        )

        np.testing.assert_array_equal(
            actual=np.array(actual['gene_list']),
            desired=np.array(expected)
        )

    else:
        with pytest.raises(err_flavor, match=err_msg):
            mapper.map_genes(
                gene_list=gene_list,
                dst_species=dst_species,
                dst_authority=dst_authority,
                ortholog_citation='NCBI',
                log=None,
                invalid_mapping_prefix=None
            )
