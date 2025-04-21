"""
Tests in this file will rely on ingestion of both our simulated
NCBI data package and our simulated bkbit data package
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
        bkbit_data_fixture0,
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

    data_file_spec = [
        {"type": "bkbit",
         "path": bkbit_data_fixture0}
    ]

    to_replace = "mmc_gene_mapper.download.download_manager.DownloadManager"
    with unittest.mock.patch(to_replace, dummy_download_mgr_fixture):
        gene_mapper = mapper.MMCGeneMapper(
            db_path=db_path,
            local_dir=tmp_dir_fixture,
            data_file_spec=data_file_spec,
            clobber=False,
            force_download=False,
            suppress_download_stdout=True
        )

    return gene_mapper


def test_get_all_authorities(mapper_fixture):
    actual = mapper_fixture.get_all_authorities()
    assert set(actual) == set(['NCBI', 'ENSEMBL'])


def test_get_all_species(mapper_fixture):
    actual = mapper_fixture.get_all_species()
    expected = [
        "mouse",
        "Mus musculus",
        "human",
        "Homo sapiens",
        "jabberwock",
        "elf",
        "fey",
        "goblin",
        "orc"
    ]
    assert set(actual) == set(expected)


def test_get_all_citations(mapper_fixture):
    actual = mapper_fixture.get_all_citations()
    assert len(actual) == 2
    expected_names = ["NCBI", "J001-2025"]
    actual_names = [citation["name"] for citation in actual]
    assert set(expected_names) == set(actual_names)


@pytest.mark.parametrize(
    "species, authority, symbol_list, expected_mapping",
    [
     ("human",
      "NCBI",
      ["symbol:0", "symbol:6", "symbol:5", "symbol:7", "nope"],
      {"symbol:0": ["NCBIGene:0"],
       "symbol:6": ["NCBIGene:6"],
       "symbol:5": [],
       "symbol:7": ["NCBIGene:5", "NCBIGene:7"],
       "nope": []
       }
      ),
     ("jabberwock",
      "ENSEMBL",
      ["symbol:24", "symbol:8", "symbol:28", "nope"],
      {"symbol:24": ["ENS26"],
       "symbol:8": [],
       "symbol:28": ["ENS30"],
       "nope": []
       }
      ),
     ("jabberwock",
      "ENSEMBL",
      ["name:24", "name:8", "name:28", "nope"],
      {"name:24": ["ENS26"],
       "name:8": [],
       "name:28": [],
       "nope": []
       }
      )
    ]
)
def test_identifiers_from_symbols(
        mapper_fixture,
        species,
        authority,
        symbol_list,
        expected_mapping):

    actual = mapper_fixture.identifiers_from_symbols(
        gene_symbol_list=symbol_list,
        species_name=species,
        authority_name=authority
    )
    assert actual['mapping'] == expected_mapping


def test_identifiers_from_symbols_error(
        mapper_fixture):

    # case where there are no citations linking a species
    # to an authority
    msg = "There are 0 citations associated with authority"
    with pytest.raises(ValueError, match=msg):
        mapper_fixture.identifiers_from_symbols(
            gene_symbol_list=["a", "b", "c"],
            species_name="human",
            authority_name="ENSEMBL"
        )

    # case where there is no such species
    msg = "no species match for flotsam"
    with pytest.raises(ValueError, match=msg):
        mapper_fixture.identifiers_from_symbols(
            gene_symbol_list=["a", "b", "c"],
            species_name="flotsam",
            authority_name="NCBI"
        )
