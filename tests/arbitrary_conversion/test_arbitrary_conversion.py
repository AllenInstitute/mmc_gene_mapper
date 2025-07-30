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
     # within human, symbol -> ENSEMBL; will not error. Could be mouse
     # but the dst_species will be used to break the tie
     (["symbol_12", "symbol_24", "symbol_74"],
      "human",
      "ENSEMBL",
      ["ENSG112", "ENSG124", "ENSG174"],
      None,
      None),
     # mouse -> mouse; use presence of mouse NCBI gene in input to
     # determine which species it is
     (["symbol_12", "symbol_24", "symbol_74", "NCBIGene:5"],
      "mouse",
      "ENSEMBL",
      ["ENSM12", "ENSM24", "ENSM74", "ENSM5"],
      None,
      None),
     # mouse -> mouse; use presence of mouse NCBI gene in input to
     # determine which species it is
     (["symbol_12", "symbol_24", "symbol_74", "NCBIGene:5"],
      "mouse",
      "NCBI",
      ["NCBIGene:12", "NCBIGene:24", "NCBIGene:74", "NCBIGene:5"],
      None,
      None),
     # mouse -> human; use presence of mouse NCBI gene in input to
     # determine which species it is
     (["symbol_12", "symbol_24", "symbol_74", "NCBIGene:5"],
      "human",
      "ENSEMBL",
      ["ENSG112", "ENSG124",
       "ortholog:UNMAPPABLE_NO_MATCH_0",
       "ortholog:UNMAPPABLE_NO_MATCH_1"],
      None,
      None),
     # mouse -> human; use presence of mouse NCBI gene in input to
     # determine which species it is
     (["symbol_12", "symbol_24", "symbol_74", "NCBIGene:5"],
      "human",
      "NCBI",
      ["NCBIGene:112", "NCBIGene:124",
       "ortholog:UNMAPPABLE_NO_MATCH_0",
       "ortholog:UNMAPPABLE_NO_MATCH_1"],
      None,
      None),
     # human -> mouse based just on symbols
     # recall can only ortholog through fish if i%2 == i%3 == 0
     # the presence of NCBIGene:136 will cue the mapper in that
     # these are human genes
     (["NCBIGene:999999", "symbol_12", "symbol_24",
       "symbol_74", "symbol_8888", "NCBIGene:136"],
      "mouse",
      "ENSEMBL",
      ["ortholog:UNMAPPABLE_NO_MATCH_0",
       "ENSM12",
       "ENSM24",
       "ortholog:UNMAPPABLE_NO_MATCH_1",
       "ortholog:UNMAPPABLE_NO_MATCH_2",
       "ENSM36"],
      None,
      None),
     # human -> mouse based just on symbols
     # recall can only ortholog through fish if i%2 == i%3 == 0
     (["NCBIGene:999999",
       "symbol_12",
       "symbol_24",
       "symbol_74",
       "symbol_6",
       "symbol_8888",
       "NCBIGene:136"],
      "mouse",
      "NCBI",
      ["ortholog:UNMAPPABLE_NO_MATCH_0",
       "NCBIGene:12",
       "NCBIGene:24",
       "ortholog:UNMAPPABLE_NO_MATCH_1",
       "NCBIGene:6",
       "ortholog:UNMAPPABLE_NO_MATCH_2",
       "NCBIGene:36"],
      None,
      None),
     # Now map from ENSEMBL to ENSEMBL
     (["NCBIGene:999999", "symbol_12", "symbol_24",
       "symbol_74", "symbol_8888", "ENSG136"],
      "mouse",
      "ENSEMBL",
      ["ortholog:UNMAPPABLE_NO_MATCH_0",
       "ENSM12",
       "ENSM24",
       "ortholog:UNMAPPABLE_NO_MATCH_1",
       "ortholog:UNMAPPABLE_NO_MATCH_2",
       "ENSM36"],
      None,
      None),
     # Now map from ENSEMBL to NCBI
     (["NCBIGene:999999",
       "symbol_12",
       "symbol_24",
       "symbol_74",
       "symbol_6",
       "symbol_8888",
       "ENSG136"],
      "mouse",
      "NCBI",
      ["ortholog:UNMAPPABLE_NO_MATCH_0",
       "NCBIGene:12",
       "NCBIGene:24",
       "ortholog:UNMAPPABLE_NO_MATCH_1",
       "NCBIGene:6",
       "ortholog:UNMAPPABLE_NO_MATCH_2",
       "NCBIGene:36"],
      None,
      None),
     # human -> mouse based just on symbols;
     # use ENSEMBL-specific symbol, which will not map to NCBI
     (["NCBIGene:999999",
       "symbol_12",
       "e_human_24",
       "symbol_74",
       "symbol_6",
       "symbol_8888",
       "NCBIGene:136"],
      "mouse",
      "NCBI",
      ["ortholog:UNMAPPABLE_NO_MATCH_0",
       "NCBIGene:12",
       "symbol:NCBI:UNMAPPABLE_NO_MATCH_0",
       "ortholog:UNMAPPABLE_NO_MATCH_2",
       "NCBIGene:6",
       "ortholog:UNMAPPABLE_NO_MATCH_3",
       "NCBIGene:36"],
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
     # -> mouse, but error because more than on species
     (["NCBIGene:112", "NCBIGene:124",
       "NCBIGene:207", "NCBIGene:205"],
      "mouse",
      "ENSEMBL",
      None,
      species_detection.InconsistentSpeciesError,
      "Multiple species inferred"
      ),
     # -> mouse, but error because more than on species
     (["NCBIGene:112", "NCBIGene:124",
       "ENSF207", "ENSF205"],
      "mouse",
      "ENSEMBL",
      None,
      species_detection.InconsistentSpeciesError,
      "ENSEMBL genes gave species"
      ),
     # again, more than one species inferred from the
     # gene identifiers that are present
     (["NCBIGene:999999", "symbol_12", "symbol_24",
       "symbol_74", "symbol_8888", "NCBIGene:136",
       "ENSM4"],
      "mouse",
      "ENSEMBL",
      None,
      species_detection.InconsistentSpeciesError,
      "ENSEMBL genes gave species"),
     # again, more than one species inferred from the
     # gene identifiers that are present
     (["NCBIGene:999999", "symbol_12", "symbol_24",
       "symbol_74", "symbol_8888", "NCBIGene:136",
       "NCBIGene:4"],
      "mouse",
      "ENSEMBL",
      None,
      species_detection.InconsistentSpeciesError,
      "Multiple species inferred"),
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
