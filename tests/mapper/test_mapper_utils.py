import pytest

import mmc_gene_mapper.mapper.mapper_utils as mapper_utils



@pytest.mark.parametrize(
    "gene_list, prefix, expected_list, expected_n",
    [(['a', 'b', 'c', 'd'], None, ['a', 'b', 'c', 'd'], 0),
     (['a', 'b', 'c', 'b', 'd', 'e'],
      None,
      ['a', 'UNMAPPABLE_DEGENERATE_0_0',
       'c', 'UNMAPPABLE_DEGENERATE_0_1',
       'd', 'e'],
      2),
     (['a', 'b', 'c', 'b', 'd', 'e'],
      'silly',
      ['a', 'silly:UNMAPPABLE_DEGENERATE_0_0',
       'c', 'silly:UNMAPPABLE_DEGENERATE_0_1',
       'd', 'e'],
      2),
     (['a', 'f', 'c', 'f', 'd', 'e', 'b', 'g', 'b'],
      None,
      ['a', 'UNMAPPABLE_DEGENERATE_1_0',
       'c', 'UNMAPPABLE_DEGENERATE_1_1',
       'd', 'e',
       'UNMAPPABLE_DEGENERATE_0_0',
       'g',
       'UNMAPPABLE_DEGENERATE_0_1'],
      4),
    ]
)
def test_mask_degenerate_ids(
        gene_list,
        prefix,
        expected_list,
        expected_n):

    (new_gene_list,
     n_degen) = mapper_utils.mask_degenerate_genes(
                     gene_list,
                     placeholder_prefix=prefix)

    assert n_degen == expected_n
    assert new_gene_list == expected_list
    assert new_gene_list is not gene_list


@pytest.mark.parametrize(
    "assign_placeholders, placeholder_prefix",
    [(True, None), (True, "test"), (False, "test")]
)
def test_apply_mapping(
        assign_placeholders,
        placeholder_prefix):

    input_genes = [
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'UNMAPPABLE_INPUT'
    ]

    mapping = {
        'a': [],
        'b': ['bb'],
        'c': ['cc1', 'cc2'],
        'd': [],
        'e': ['ee'],
        'f': ['ff'],
        'g': ['ee'],
        'h': ['hh'],
        'i': ['ff']
    }
    result = mapper_utils.apply_mapping(
        gene_list=input_genes,
        mapping=mapping,
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )
    expected_failure_log = {
            'zero matches': 3,
            'many matches': 1,
            'degenerate matches': 4
        }

    if assign_placeholders:
        if placeholder_prefix is None:
            expected_gene_list = [
                'UNMAPPABLE_NO_MATCH_0',
                'bb',
                'UNMAPPABLE_MANY_MATCHES_0',
                'UNMAPPABLE_NO_MATCH_1',
                'UNMAPPABLE_DEGENERATE_0_0',
                'UNMAPPABLE_DEGENERATE_1_0',
                'UNMAPPABLE_DEGENERATE_0_1',
                'hh',
                'UNMAPPABLE_DEGENERATE_1_1',
                'UNMAPPABLE_INPUT'
            ]
        else:
            expected_gene_list = [
                'test:UNMAPPABLE_NO_MATCH_0',
                'bb',
                'test:UNMAPPABLE_MANY_MATCHES_0',
                'test:UNMAPPABLE_NO_MATCH_1',
                'test:UNMAPPABLE_DEGENERATE_0_0',
                'test:UNMAPPABLE_DEGENERATE_1_0',
                'test:UNMAPPABLE_DEGENERATE_0_1',
                'hh',
                'test:UNMAPPABLE_DEGENERATE_1_1',
                'UNMAPPABLE_INPUT'
            ]
    else:
        expected_gene_list = [
            'a', 'bb', 'c', 'd',
            'test:UNMAPPABLE_DEGENERATE_0_0',
            'test:UNMAPPABLE_DEGENERATE_1_0',
            'test:UNMAPPABLE_DEGENERATE_0_1',
            'hh',
            'test:UNMAPPABLE_DEGENERATE_1_1',
            'UNMAPPABLE_INPUT'
        ]

    expected = {
        'failure_log': expected_failure_log,
        'gene_list': expected_gene_list
    }

    assert result == expected


def test_apply_mapping_override_degenerate_placeholder():
    """
    Test that two genes which come in with degenerate placeholder
    identifiers get mapped to unique identifiers
    """
    gene_list = ['a', 'UNMAPPABLE', 'b']
    mapping = {
        'a': ['aa'],
        'b': ['bb']
    }
    actual = mapper_utils.apply_mapping(
        gene_list=gene_list,
        mapping=mapping,
        assign_placeholders=False,
        placeholder_prefix=None
    )
    assert actual['gene_list'] == ['aa', 'UNMAPPABLE', 'bb']

    gene_list = ['a', 'UNMAPPABLE', 'b', 'UNMAPPABLE']
    actual = mapper_utils.apply_mapping(
        gene_list=gene_list,
        mapping=mapping,
        assign_placeholders=False,
        placeholder_prefix=None
    )
    assert actual['gene_list'] == [
        'aa',
        'UNMAPPABLE_DEGENERATE_0_0',
        'bb',
        'UNMAPPABLE_DEGENERATE_0_1'
    ]

    gene_list = ['a', 'UNMAPPABLE', 'b', 'UNMAPPABLE']
    actual = mapper_utils.apply_mapping(
        gene_list=gene_list,
        mapping=mapping,
        assign_placeholders=False,
        placeholder_prefix='silly'
    )
    assert actual['gene_list'] == [
        'aa',
        'silly:UNMAPPABLE_DEGENERATE_0_0',
        'bb',
        'silly:UNMAPPABLE_DEGENERATE_0_1'
    ]
