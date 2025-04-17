"""
Ingest a mock NCBI package and try to query it for
equivalences and orthologs
"""

import pytest

import tempfile
import unittest.mock

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.mapper.mapper as mapper


@pytest.fixture(scope='session')
def mapper_fixture(
        ncbi_file_package_fixture,
        dummy_download_mgr_fixture,
        tmp_dir_fixture):
    """
    Return an instantiation of the MMCGeneMapper class
    based on our simulated NCBI file package
    """
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_fixture)
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir,
        prefix='ncbi_gene_mapper_',
        suffix='.db',
        delete=True
    )
    to_replace = "mmc_gene_mapper.download.download_manager.DownloadManager"
    with unittest.mock.patch(to_replace, dummy_download_mgr_fixture):
        gene_mapper = mapper.MMCGeneMapper(
            db_path=db_path,
            local_dir=tmp_dir_fixture,
            data_file_spec=None,
            clobber=False,
            force_download=False,
            suppress_download_stdout=True
        )

    return gene_mapper


def test_get_orthologs_from_ncbi(
        mapper_fixture):

    gene_idx_list = [0, 1, 2, 4, 7, 6]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapper_fixture.ortholog_genes(
        authority='NCBI',
        src_species_name='human',
        dst_species_name='jabberwock',
        gene_list=gene_list,
        citation_name='NCBI'
    )

    expected = {
        'NCBIGene:0': ['NCBIGene:14'],
        'NCBIGene:1': ['NCBIGene:13'],
        'NCBIGene:2': [],
        'NCBIGene:4': ['NCBIGene:12'],
        'NCBIGene:6': [],
        'NCBIGene:7': ['NCBIGene:15']
    }
    assert actual['mapping'] == expected

    gene_idx_list = [20, 21, 22, 23, 24, 27]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapper_fixture.ortholog_genes(
        authority='NCBI',
        src_species_name='mouse',
        dst_species_name='jabberwock',
        gene_list=gene_list,
        citation_name='NCBI'
    )

    expected = {
        'NCBIGene:20': [],
        'NCBIGene:21': ['NCBIGene:14'],
        'NCBIGene:22': [],
        'NCBIGene:23': ['NCBIGene:12'],
        'NCBIGene:24': [],
        'NCBIGene:27': ['NCBIGene:13']
    }
    assert actual['mapping'] == expected


@pytest.mark.skip('cannot get equivalent genes until we ingest ENSEMBL data')
def test_get_equivalent_genes_from_ncbi(
        mapper_fixture):

    gene_idx_list = [1, 2, 3, 5, 6, 7, 8]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapper_fixture.equivalent_genes(
        input_authority='NCBI',
        output_authority='ENSEMBL',
        gene_list=gene_list,
        species_name='human',
        citation_name='NCBI'
    )
