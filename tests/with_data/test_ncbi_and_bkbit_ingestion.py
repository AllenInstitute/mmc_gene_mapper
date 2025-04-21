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
        alt_ortholog_file_fixture,
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
         "path": bkbit_data_fixture0
         },
        {"type": "hmba_orthologs",
         "path": alt_ortholog_file_fixture,
         "name": "alternative_orthologs"
         }
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
    assert len(actual) == 3
    expected_names = ["NCBI", "J001-2025", "alternative_orthologs"]
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

    gene_idx_list = [20, 21, 22, 23, 24, 27]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapper_fixture.ortholog_genes(
        authority='NCBI',
        src_species_name='mouse',
        dst_species_name='human',
        gene_list=gene_list,
        citation_name='NCBI'
    )

    expected = {
        'NCBIGene:20': [],
        'NCBIGene:21': ['NCBIGene:0'],
        'NCBIGene:22': [],
        'NCBIGene:23': ['NCBIGene:4'],
        'NCBIGene:24': ['NCBIGene:6'],
        'NCBIGene:27': ['NCBIGene:1']
    }
    assert actual['mapping'] == expected


def test_get_equivalent_genes_from_ncbi(
        mapper_fixture):

    gene_idx_list = [1, 2, 3, 6, 14, 10]
    gene_list = [f'ENS{ii}' for ii in gene_idx_list]
    actual = mapper_fixture.equivalent_genes(
        input_authority='ENSEMBL',
        output_authority='NCBI',
        gene_list=gene_list,
        species_name='human',
        citation_name='NCBI'
    )

    expected = {
        'ENS1': [],
        'ENS2': ['NCBIGene:1'],
        'ENS3': [],
        'ENS6': ['NCBIGene:3'],
        'ENS14': ['NCBIGene:5', 'NCBIGene:7'],
        'ENS10': ['NCBIGene:5']
    }

    assert set(expected.keys()) == set(actual['mapping'].keys())
    for k in expected:
        assert set(expected[k]) == set(actual['mapping'][k])


@pytest.mark.parametrize(
    "src_species, dst_species, citation, gene_list, expected_mapping",
    [("human",
      "jabberwock",
      "alternative_orthologs",
      ["NCBIGene:0", "NCBIGene:2", "NCBIGene:7", "NCBIGene:4", "NCBIGene:6"],
      {"NCBIGene:0": ["NCBIGene:11"],
       "NCBIGene:2": ["NCBIGene:13"],
       "NCBIGene:7": ["NCBIGene:12"],
       "NCBIGene:4": [],
       "NCBIGene:6": []
       }
      ),
     ("mouse",
      "human",
      "alternative_orthologs",
      ["NCBIGene:21", "NCBIGene:22", "NCBIGene:23",
       "NCBIGene:24", "NCBIGene:25", "NCBIGene:26",
       "NCBIGene:27"],
      {"NCBIGene:21": [],
       "NCBIGene:22": [],
       "NCBIGene:23": [],
       "NCBIGene:24": [],
       "NCBIGene:25": ["NCBIGene:0"],
       "NCBIGene:26": ["NCBIGene:2"],
       "NCBIGene:27": []
       }
      ),
     ("mouse",
      "human",
      "NCBI",
      ["NCBIGene:21", "NCBIGene:22", "NCBIGene:23",
       "NCBIGene:24", "NCBIGene:25", "NCBIGene:26",
       "NCBIGene:27"],
      {"NCBIGene:21": ["NCBIGene:0"],
       "NCBIGene:22": [],
       "NCBIGene:23": ["NCBIGene:4"],
       "NCBIGene:24": ["NCBIGene:6"],
       "NCBIGene:25": [],
       "NCBIGene:26": [],
       "NCBIGene:27": ["NCBIGene:1"]
       }
      ),
     ("mouse",
      "jabberwock",
      "NCBI",
      ["NCBIGene:21", "NCBIGene:22", "NCBIGene:23",
       "NCBIGene:24", "NCBIGene:25", "NCBIGene:26",
       "NCBIGene:27"],
      {"NCBIGene:21": ["NCBIGene:14"],
       "NCBIGene:22": [],
       "NCBIGene:23": ["NCBIGene:12"],
       "NCBIGene:24": [],
       "NCBIGene:25": [],
       "NCBIGene:26": [],
       "NCBIGene:27": ["NCBIGene:13"]
       }
      ),
     ]
)
def test_alternative_orthologs(
        mapper_fixture,
        src_species,
        dst_species,
        citation,
        gene_list,
        expected_mapping):

    actual = mapper_fixture.ortholog_genes(
        authority="NCBI",
        src_species_name=src_species,
        dst_species_name=dst_species,
        citation_name=citation,
        gene_list=gene_list
    )
    assert actual['mapping'] == expected_mapping
