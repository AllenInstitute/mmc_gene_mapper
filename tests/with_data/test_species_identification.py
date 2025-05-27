import pytest

import numpy as np

import mmc_gene_mapper.mapper.mapper_utils as mapper_utils


def test_determine_species_and_authority_from_data(
        mapper_fixture):

    garbage_genes = ["a", "b", "c", "d"]
    actual = mapper_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=garbage_genes
    )

    assert actual['species'] is None
    assert actual['species_taxon'] is None
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['symbol']*len(garbage_genes))
    )

    jabberwock_genes = ["ENSX22", "ENSX26"]
    actual = mapper_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    print(actual)
    assert actual['species'] == 'jabberwock'
    assert actual['species_taxon'] == 999
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['ENSEMBL', 'ENSEMBL'])
    )

    jabberwock_genes = ["NCBIGene:11", "NCBIGene:13"]
    actual = mapper_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    assert actual['species'] == 'jabberwock'
    assert actual['species_taxon'] == 999
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['NCBI', 'NCBI'])
    )

    jabberwock_genes = ["a0", "a1", "NCBIGene:11", "a2",
                        "ENSX22", "NCBIGene:13", "ENSX26"]
    actual = mapper_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )

    assert actual['species'] == 'jabberwock'
    assert actual['species_taxon'] == 999
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(
            ['symbol',
             'symbol',
             'NCBI',
             'symbol',
             'ENSEMBL',
             'NCBI',
             'ENSEMBL'
             ]
        )
    )

    mixed_genes = ["a2", "NCBIGene:11", "a3", "ENSX2", "NCBIGene:13", "a4"]
    msg = (
        "NCBI genes gave species 'jabberwock'"
    )
    with pytest.raises(mapper_utils.InconsistentSpeciesError, match=msg):
        mapper_utils.detect_species_and_authority(
            db_path=mapper_fixture.db_path,
            gene_list=mixed_genes
        )


def test_backend_determine_species_and_authority_from_data(
        mapper_fixture):

    garbage_genes = ["a", "b", "c", "d"]
    actual = mapper_utils._detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=garbage_genes
    )
    expected = {
        "authority": None,
        "species": None,
        "species_taxon": None}
    assert actual == expected

    jabberwock_genes = ["ENSX22", "ENSX26"]
    actual = mapper_utils._detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    expected = {
        "authority": "ENSEMBL",
        "species": "jabberwock",
        "species_taxon": 999}
    assert actual == expected

    jabberwock_genes = [f'a{ii}' for ii in range(45)] + ["ENSX22", "ENSX26"]
    assert len(jabberwock_genes) > 25
    actual = mapper_utils._detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    expected = {
        "authority": "ENSEMBL",
        "species": "jabberwock",
        "species_taxon": 999}
    assert actual == expected

    jabberwock_genes = ["NCBIGene:11", "NCBIGene:13"]
    actual = mapper_utils._detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    expected = {
        "authority": "NCBI",
        "species": "jabberwock",
        "species_taxon": 999
    }
    assert actual == expected

    msg = "Multiple authorities inferred"
    mixed_genes = ["NCBIGene:11", "ENSX22"]
    with pytest.raises(mapper_utils.InconsistentSpeciesError, match=msg):
        mapper_utils._detect_species_and_authority(
            db_path=mapper_fixture.db_path,
            gene_list=mixed_genes
        )

    msg = "Multiple species inferred"
    mixed_genes = ["NCBIGene:11", "ENSX22", "NCBIGene:3"]
    with pytest.raises(mapper_utils.InconsistentSpeciesError, match=msg):
        mapper_utils._detect_species_and_authority(
            db_path=mapper_fixture.db_path,
            gene_list=mixed_genes
        )
