import pytest

import numpy as np

import mmc_gene_mapper.utils.str_utils as str_utils


def test_int_from_identifier():
    assert str_utils.int_from_identifier('ENS009991') == 9991
    assert str_utils.int_from_identifier('ENSX17') == 17
    assert str_utils.int_from_identifier('ABC89910') == 89910
    assert str_utils.int_from_identifier('AB123DEF8910') == 8910

    with pytest.raises(ValueError, match='could not get one integer'):
        str_utils.int_from_identifier('abcabc')

    with pytest.raises(ValueError, match='could not get one integer'):
        str_utils.int_from_identifier('abc123def')


def test_characterize_gene_identifiers():

    values = [
       'a',
       '112321',
       'NCBI:5513',
       'ENSMUSG00998811',
       'ENST778881123',
       'ENST78132b',
       'NCBIGene:7781',
       'NCBIGene::7777123',
       'ENSMUSG',
       'NCBI',
       'abcNCBIGene:777',
       'abcENSX66'
    ]
    expected = [
        'symbol',
        'symbol',
        'NCBI',
        'ENSEMBL',
        'ENSEMBL',
        'symbol',
        'NCBI',
        'symbol',
        'symbol',
        'symbol',
        'symbol',
        'symbol'
    ]

    actual = str_utils.characterize_gene_identifiers(
        values
    )

    np.testing.assert_array_equal(
        actual,
        expected
    )
