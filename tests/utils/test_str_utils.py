import mmc_gene_mapper.utils.str_utils as str_utils


def test_detect_gene_identifiers():

    assert str_utils.detect_gene_identifiers(['a', 'b', 'c']) is None
    assert str_utils.detect_gene_identifiers(
        ['NCBIGene:011', 'NCBIGene:9191', 'NCBIGene:223']) == 'NCBI'
    assert str_utils.detect_gene_identifiers(
        ['ENSMUSG0001', 'ENSG88881', 'ENST888111']) == 'ENSEMBL'
    assert str_utils.detect_gene_identifiers(
        ['ENSMUST0091', 'NCBIGene:99']) is None
