import pytest

import mmc_gene_mapper.mapper.mapper_utils as mapper_utils


@pytest.mark.parametrize(
    "assign_placeholders, placeholder_prefix",
    [(True, None), (True, "test"), (False, "test")]
)
def test_apply_mapping(
        assign_placeholders,
        placeholder_prefix):

    input_genes = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
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
            'zero matches': 2,
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
                'UNMAPPABLE_DEGENERATE_1_1'
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
                'test:UNMAPPABLE_DEGENERATE_1_1'
            ]
    else:
        expected_gene_list = [
            'a', 'bb', 'c', 'd',
            'test:UNMAPPABLE_DEGENERATE_0_0',
            'test:UNMAPPABLE_DEGENERATE_1_0',
            'test:UNMAPPABLE_DEGENERATE_0_1',
            'hh',
            'test:UNMAPPABLE_DEGENERATE_1_1'
        ]

    expected = {
        'failure_log': expected_failure_log,
        'gene_list': expected_gene_list
    }

    assert result == expected
