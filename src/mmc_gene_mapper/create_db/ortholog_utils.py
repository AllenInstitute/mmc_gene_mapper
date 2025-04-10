import numpy as np

def assign_ortholog_group(gene0_list, gene1_list):
    filtered0 = []
    filtered1 = []
    for g0, g1 in zip(gene0_list, gene1_list):
        if g0 != g1:
            filtered0.append(g0)
            filtered1.append(g1)

    graph = create_ortholog_graph(
        gene0_list=filtered0,
        gene1_list=filtered1)
    print("    GOT GRAPH")
    return assign_ortholog_group_from_graph(
        graph=graph,
        root_gene_list=sorted(set(filtered0)))


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


def assign_ortholog_group_from_graph(
        graph,
        root_gene_list=None):
    """
    Assign genes to ortholog groups (a group is
    a numerical index used to indicate genes that are
    all orthologs of each other)

    Parameters
    ----------
    graph:
        a graph as produced by ortholog_utils.create_ortholog_graph
    root_gene_list:
        a list of genes to start walking from first (because we think
        they will produce the largest initial groups)

    Returns
    -------
    A dict mapping each gene in the graph to its
    ortholog group.

    Note
    ----
    graph gets slowly whittled down in-place to save time
    searching over nodes.
    """
    if root_gene_list is None:
        root_gene_list = []
    else:
        root_gene_list = np.array(root_gene_list)
        n_neigh = np.array([len(graph[g]) for g in root_gene_list])
        sorted_dex = np.argsort(n_neigh)[-1::-1]
        root_gene_list = list(root_gene_list[sorted_dex])

    group_idx = 0
    gene_to_group = dict()
    n0 = len(graph)
    last_decrement = 0
    while True:
        start_node = None
        for ii, node in enumerate(root_gene_list):
            if node in graph:
                start_node = node
                to_pop = ii
                break
        if start_node is not None:
            root_gene_list.pop(to_pop)
        else:
            if len(graph) == 0:
                break
            start_node = list(graph.keys())[0]

        (these_nodes) = _walk_from(
             graph=graph,
             start_node=start_node)

        for node in these_nodes:
            assert node not in gene_to_group
            gene_to_group[node] = group_idx

        group_idx += 1
        n1 = len(graph)
        decrement = n0-n1
        if group_idx == 1 or decrement > last_decrement*10:
            print(f"    group {group_idx} decrement {decrement:.3e} of {n0:.3e}")
            last_decrement = decrement

    return gene_to_group


def _walk_from(graph, start_node):
    """
    Walk a graph starting from start_node. Return the set
    of all nodes connected to start_node
    """

    queue = []
    already_visited = set()
    already_visited.add(start_node)
    for neigh in graph.pop(start_node):
        queue.append(neigh)

    while len(queue) > 0:
        chosen_idx = queue[0]
        queue = queue[1:]
        if chosen_idx in already_visited:
            continue
        chosen_neigh = graph.pop(chosen_idx)
        chosen_neigh = list(chosen_neigh-already_visited)
        already_visited.add(chosen_idx)
        queue += chosen_neigh

    return already_visited
