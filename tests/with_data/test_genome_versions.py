"""
Tests in this module will cover the case where, for instance,
the NCBI data references ENSEMBL genes that are not in the ENSEMBL data
"""
import pytest

import numpy as np


@pytest.mark.parametrize(
    'gene_list, expected',
    [
     (['NCBIGene:24', 'NCBIGene:25'],
      ['NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_0',
       'ENSX50'
       ]
      ),
     (['NCBIGene:24', 'NCBIGene:25', 'NCBIGene:31'],
      ['NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_0',
       'ENSX50',
       'NCBI:ENSEMBL:UNMAPPABLE_NO_MATCH_1'
       ]
      ),
    ]
)
def test_mapping_to_when_ensembl_missing(
        mapper_fixture,
        gene_list,
        expected):
    result = mapper_fixture.map_genes(
        gene_list=gene_list,
        dst_authority='ENSEMBL',
        dst_species='Mus musculus'
    )

    np.testing.assert_array_equal(
        np.array(result['gene_list']),
        np.array(expected)
    )
