"""
Test that the metadata returned with a mapping records
all of the necessary steps
"""
import pytest

import sqlite3



def test_mapper_metadata(
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
