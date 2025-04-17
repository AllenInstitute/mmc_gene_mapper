def get_gene_to_species_map(
        cursor,
        gene_list,
        authority_idx,
        chunk_size=10000):
    """
    Return a dict mapping gene IDs to species IDs
    according to a specified authority.

    Parameters
    ----------
    gene_list:
        list of integers specifying the genes
    authority_idx:
        integer specifying the authority
    chunk_size:
        number of genes to query at a time

    Returns
    -------
    A dict mapping gene ID to species ID (both integers)

    Notes
    -----
    The function will return all (authority_idx, citation) combinations
    for a given species. If the gene is assigned to different species
    accross citations (which should not be possible in well-formed
    data) an error will be thrown.

    Any genes not listed in the gene table simply will not appear
    in the gene_to_species map.
    """
    gene_to_species_taxon = dict()
    for i0 in range(0, len(gene_list), chunk_size):
        gene_subset = gene_list[i0:i0+chunk_size]
        query = """
            SELECT
                id,
                species_taxon
            FROM gene
            WHERE
                authority=?
            AND
                id IN (
        """
        query += ",".join(["?"]*len(gene_subset))
        query += ")"
        raw = cursor.execute(
            query,
            (authority_idx,
             *gene_subset)
        ).fetchall()
        for row in raw:
            if row[0] in gene_to_species_taxon:
                if gene_to_species_taxon[row[0]] != row[1]:
                    raise ValueError(
                        "Conflicting species taxon for "
                        f"gene_id {row[0]} authority {authority}; "
                        "unclear how to proceed."
                    )
            gene_to_species_taxon[row[0]] = row[1]
    return gene_to_species_taxon
