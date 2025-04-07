import numpy as np

import mmc_gene_mapper.utils.str_utils as str_utils


def test_characterize_gene_identifiers():

    values = [
       'a',
       '112321',
       'NCBI:5513',
       'ENSMUSG00998811',
       'ENST778881123',
       'ENST78132b',
       'NCBIGene:7781',
       'NCBIGene::7777123'
    ]
    expected = [
        'symbol',
        'symbol',
        'NCBI',
        'ENSEMBL',
        'ENSEMBL',
        'symbol',
        'NCBI',
        'symbol'
    ]

    actual = str_utils.characterize_gene_identifiers(
        values
    )

    np.testing.assert_array_equal(
        actual,
        expected
    )
