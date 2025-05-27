import pytest

import mmc_gene_mapper.mapper.mapper_utils as mapper_utils


def test_determine_species_and_authority_from_data(
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

    jabberwock_genes = ["ENS22", "ENS26"]
    actual = mapper_utils._detect_species_and_authority(
        db_path=mapper_fixture.db_path,
        gene_list=jabberwock_genes
    )
    expected = {
        "authority": "ENSEMBL",
        "species": "jabberwock",
        "species_taxon": 999}
    assert actual == expected

    jabberwock_genes = [f'a{ii}' for ii in range(45)] + ["ENS22", "ENS26"]
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
    mixed_genes = ["NCBIGene:11", "ENS22"]
    with pytest.raises(mapper_utils.InconsistentSpeciesError, match=msg):
        mapper_utils._detect_species_and_authority(
            db_path=mapper_fixture.db_path,
            gene_list=mixed_genes
        )

    msg = "Multiple species inferred"
    mixed_genes = ["NCBIGene:11", "ENS22", "NCBIGene:3"]
    with pytest.raises(mapper_utils.InconsistentSpeciesError, match=msg):
        mapper_utils._detect_species_and_authority(
            db_path=mapper_fixture.db_path,
            gene_list=mixed_genes
        )
