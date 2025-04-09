import pytest

import mmc_gene_mapper.create_db.ortholog_utils as ortholog_utils


def test_create_ortholog_graph():

    gene0_list = [
        1, 1, 1, 2, 2, 3, 3, 4
    ]
    gene1_list = [
        2, 3, 5, 4, 6, 7, 8, 5
    ]

    graph = ortholog_utils.create_ortholog_graph(
        gene0_list=gene0_list,
        gene1_list=gene1_list
    )
    assert len(graph) == 8
    assert set(graph.keys()) == set(range(1, 9, 1))

    assert graph[1] == set([2, 3, 5])
    assert graph[2] == set([1, 4, 6])
    assert graph[3] == set([1, 7, 8])
    assert graph[4] == set([2, 5])
    assert graph[5] == set([1, 4])
    assert graph[6] == set([2])
    assert graph[7] == set([3])
    assert graph[8] == set([3])

    with pytest.raises(ValueError, match="must be the same size"):
        ortholog_utils.create_ortholog_graph(
            gene0_list=[1, 2, 3],
            gene1_list=[4, 5]
        )
