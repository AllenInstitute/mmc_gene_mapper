import mmc_gene_mapper.mapper.mapper_utils as mapper_utils


def test_apply_mapping():

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
        mapping=mapping
    )

    expected = {
        'failure_log': {
            'zero matches': 2,
            'many matches': 1,
            'degenerate matches': 4
        },
        'gene_list': [
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
    }

    assert result == expected
