import pytest

import copy
import mmc_gene_mapper.create_db.ortholog_utils as ortholog_utils


def test_create_ortholog_graph():

    gene0_list = [
        1, 1, 1, 2, 2, 3, 3, 4, 9
    ]
    gene1_list = [
        2, 3, 5, 4, 6, 7, 8, 5, 10
    ]

    graph = ortholog_utils.create_ortholog_graph(
        gene0_list=gene0_list,
        gene1_list=gene1_list
    )
    assert len(graph) == 10
    assert set(graph.keys()) == set(range(1, 11, 1))

    assert graph[1] == set([2, 3, 5])
    assert graph[2] == set([1, 4, 6])
    assert graph[3] == set([1, 7, 8])
    assert graph[4] == set([2, 5])
    assert graph[5] == set([1, 4])
    assert graph[6] == set([2])
    assert graph[7] == set([3])
    assert graph[8] == set([3])
    assert graph[9] == set([10])
    assert graph[10] == set([9])

    with pytest.raises(ValueError, match="must be the same size"):
        ortholog_utils.create_ortholog_graph(
            gene0_list=[1, 2, 3],
            gene1_list=[4, 5]
        )


@pytest.mark.parametrize("use_guess", [True, False])
def test_assign_ortholog_group_from_graph(use_guess):

    graph = dict()
    graph[1] = set([2, 3, 5])
    graph[2] = set([1, 4, 6])
    graph[3] = set([1, 7, 8])
    graph[4] = set([2, 5])
    graph[5] = set([1, 4])
    graph[6] = set([2])
    graph[7] = set([3])
    graph[8] = set([3])
    graph[9] = set([10])
    graph[10] = set([9])

    if use_guess:
        root_gene_list = [1, 2, 5]
    else:
        root_gene_list = None

    group_lookup = ortholog_utils.assign_ortholog_group_from_graph(
        graph=copy.deepcopy(graph),
        root_gene_list=root_gene_list)

    for k in graph:
        assert k in group_lookup

    for k in [1, 2, 3, 4, 5, 6, 7, 8]:
        assert group_lookup[k] == group_lookup[1]
    assert group_lookup[9] == group_lookup[10]
    assert group_lookup[9] != group_lookup[1]


def test_assign_ortholog_group():
    gene0_list = [
        1, 1, 1, 2, 2, 3, 3, 4, 9
    ]
    gene1_list = [
        2, 3, 5, 4, 6, 7, 8, 5, 10
    ]

    group_lookup = ortholog_utils.assign_ortholog_group(
        gene0_list=gene0_list,
        gene1_list=gene1_list)

    for k in range(1, 11,1):
        assert k in group_lookup

    for k in [1, 2, 3, 4, 5, 6, 7, 8]:
        assert group_lookup[k] == group_lookup[1]
    assert group_lookup[9] == group_lookup[10]
    assert group_lookup[9] != group_lookup[1]
