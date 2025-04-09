

def create_ortholog_graph(gene0_list, gene1_list):
    """
    Take two lists of genes such that
    gene0_list[ii] is an ortholog of gene1_list[ii]
    construct a graph that links all ortholog groups.

    The graph will be a dict. Keys are genes. Values are
    sets of genes that are linked directly to the keys in
    the graph. Each (gene0, gene1) pair is entered twice
    so that the graph is bi-directional
    """
    if len(gene0_list) != len(gene1_list):
        raise ValueError(
            f"gene0_list has {len(gene0_list)} elements; "
            f"gene1_list has {len(gene1_list)} elements; "
            "these must be the same size"
        )
    graph = dict()
    for g0, g1 in zip(gene0_list, gene1_list):
        g0 = int(g0)
        g1 = int(g1)
        if g0 not in graph:
            graph[g0] = set()
        graph[g0].add(g1)
        if g1 not in graph:
            graph[g1] = set()
        graph[g1].add(g0)

    return graph
