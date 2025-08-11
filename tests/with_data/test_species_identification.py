import numpy as np

import mmc_gene_mapper.metadata.classes as metadata_classes
import mmc_gene_mapper.mapper.species_detection as species_utils


def test_determine_species_and_authority_from_data(
        mapper_fixture):

    garbage_genes = ["a", "b", "c", "d"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=garbage_genes
    )

    assert actual['species'] is None
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['symbol']*len(garbage_genes))
    )

    jabberwock_genes = ["ENSX22", "ENSX26"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    assert isinstance(actual['species'], metadata_classes.Species)
    assert actual['species'].name == 'jabberwock'
    assert actual['species'].taxon == 999
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['ENSEMBL', 'ENSEMBL'])
    )

    jabberwock_genes = ["a0", "ENSX22", "a1", "ENSX26"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    assert actual['species'].name == 'jabberwock'
    assert actual['species'].taxon == 999
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['symbol', 'ENSEMBL', 'symbol', 'ENSEMBL'])
    )

    # test case where ENS gives a species and NCBI gives None
    jabberwock_genes = ["a0", "ENSX22", "a1", "ENSX26", "NCBIGene:1111111"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    assert actual['species'].name == 'jabberwock'
    assert actual['species'].taxon == 999
    np.testing.assert_array_equal(
        actual=actual['authority'],
        desired=np.array(['symbol', 'ENSEMBL', 'symbol', 'ENSEMBL', 'symbol'])
    )

    jabberwock_genes = ["NCBIGene:11", "NCBIGene:13"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    assert actual['species'].name == 'jabberwock'
    assert actual['species'].taxon == 999
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['NCBI', 'NCBI'])
    )

    # test case where NCBI gives non-None and ENSEMBL gives None
    jabberwock_genes = ["NCBIGene:11", "NCBIGene:13", "ENSX777777777"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    assert actual['species'].name == 'jabberwock'
    assert actual['species'].taxon == 999
    np.testing.assert_array_equal(
        actual['authority'],
        np.array(['NCBI', 'NCBI', 'symbol'])
    )

    jabberwock_genes = ["a0", "a1", "NCBIGene:11", "a2",
                        "ENSX22", "NCBIGene:13", "ENSX26"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )

    assert actual['species'].name == 'jabberwock'
    assert actual['species'].taxon == 999
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

    # voting should choose jabberwock
    mixed_genes = ["a2", "NCBIGene:11", "a3", "ENSX2", "NCBIGene:13", "a4"]
    actual = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=mixed_genes
    )
    assert actual['species'].name == 'jabberwock'


def test_species_from_odd_symbols(mapper_fixture):
    """
    Test that symbols which fit the ENSEMBL and NCBI regex
    return species==None
    """
    ens_list = [
        'ENSBA0',
        'ENSCA1'
    ]
    result = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=ens_list
    )
    assert result['species'] is None
    np.testing.assert_array_equal(
        result['authority'],
        np.array(['symbol']*2)
    )

    ncbi_list = [
        'NCBIBA0',
        'NCBICA1'
    ]
    result = species_utils.detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=ncbi_list
    )
    assert result['species'] is None
    np.testing.assert_array_equal(
        result['authority'],
        np.array(['symbol']*2)
    )
