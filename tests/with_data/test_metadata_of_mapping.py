"""
Test that the metadata returned with a mapping records
all of the necessary steps
"""


def test_mapper_metadata_full_suite(
        mapper_fixture):

    gene_list = [
        'symbol:0',
        'NCBIGene:1',
        'ENSX6',
        'symbol:2',
        'symbol:4',
        'NCBIGene:5',
        'NCBIGene:6',
        'NCBIGene:7',
        'symbol:8',
        'ENSX18'
    ]

    result = mapper_fixture.map_genes(
        gene_list=gene_list,
        dst_species='jabberwock',
        dst_authority='ENSEMBL',
        ortholog_citation='NCBI')

    assert len(result) == 2
    assert isinstance(result, dict)
    assert 'metadata' in result
    assert 'gene_list' in result

    metadata = result['metadata']
    assert len(metadata) == 4

    # check metadata for initial mapping to NCBI
    first_step = metadata[:2]
    found_ens_to_ncbi = False
    found_sym_to_ncbi = False
    for step in first_step:
        assert 'citation' in step
        mapping = step['mapping']
        assert mapping['axis'] == 'authority'
        assert mapping['to'] == 'NCBI'
        if mapping['from'] == 'symbol':
            found_sym_to_ncbi = True
        elif mapping['from'] == 'ENSEMBL':
            found_ens_to_ncbi = True
    assert found_sym_to_ncbi
    assert found_ens_to_ncbi

    # check mapping from species to species
    assert 'citation' in metadata[2]
    mapping = metadata[2]['mapping']
    assert mapping['axis'] == 'species'
    assert mapping['from'] == {'name': 'human', 'taxon': 9606}
    assert mapping['to'] == {'name': 'jabberwock', 'taxon': 999}

    # check final mapping to ENSEMBL
    assert 'citation' in metadata[3]
    mapping = metadata[3]['mapping']
    assert mapping['axis'] == 'authority'
    assert mapping['from'] == 'NCBI'
    assert mapping['to'] == 'ENSEMBL'


def test_mapper_metadata_no_ortholog(
        mapper_fixture):

    gene_list = [
        'NCBIGene:0',
        'NCBIGene:1',
        'ENSX6',
        'NCBIGene:2',
        'NCBIGene:4',
        'NCBIGene:5',
        'NCBIGene:6',
        'NCBIGene:7',
        'NCBIGene:8',
        'ENSX18'
    ]

    result = mapper_fixture.map_genes(
        gene_list=gene_list,
        dst_species='human',
        dst_authority='ENSEMBL',
        ortholog_citation='NCBI')

    assert len(result) == 2
    assert isinstance(result, dict)
    assert 'metadata' in result
    assert 'gene_list' in result

    metadata = result['metadata']
    assert len(metadata) == 1
    assert 'citation' in metadata[0]
    mapping = metadata[0]['mapping']
    assert mapping['axis'] == 'authority'
    assert mapping['to'] == 'ENSEMBL'
    assert mapping['from'] == 'NCBI'


def test_mapper_metadata_just_ortholog(
        mapper_fixture):

    gene_list = [
        'NCBIGene:0',
        'NCBIGene:1',
        'NCBIGene:3',
        'NCBIGene:2',
        'NCBIGene:4',
        'NCBIGene:5',
        'NCBIGene:6',
        'NCBIGene:7',
        'NCBIGene:8',
        'NCBIGene:9'
    ]

    result = mapper_fixture.map_genes(
        gene_list=gene_list,
        dst_species='mouse',
        dst_authority='NCBI',
        ortholog_citation='NCBI')

    assert len(result) == 2
    assert isinstance(result, dict)
    assert 'metadata' in result
    assert 'gene_list' in result

    metadata = result['metadata']
    assert len(metadata) == 1
    assert 'citation' in metadata[0]
    mapping = metadata[0]['mapping']
    assert mapping['axis'] == 'species'
    assert mapping['to'] == {'name': 'mouse', 'taxon': 10090}
    assert mapping['from'] == {'name': 'human', 'taxon': 9606}
