"""
Tests in this file will rely on ingestion of both our simulated
NCBI data package and our simulated bkbit data package
"""
import pytest

import sqlite3

import mmc_gene_mapper.query_db.query as query_utils
import mmc_gene_mapper.mapper.mapping_functions as mapping_functions


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
    assert len(actual) == 4
    expected_names = [
        "NCBI",
        "J001-2025",
        "M001-2025",
        "alternative_orthologs"]
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
     ("9606",
      "NCBI",
      ["symbol:0", "symbol:6", "symbol:5", "symbol:7", "nope"],
      {"symbol:0": ["NCBIGene:0"],
       "symbol:6": ["NCBIGene:6"],
       "symbol:5": [],
       "symbol:7": ["NCBIGene:5", "NCBIGene:7"],
       "nope": []
       }
      ),
     (9606,
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
      {"symbol:24": ["ENSX26"],
       "symbol:8": [],
       "symbol:28": ["ENSX30"],
       "nope": []
       }
      ),
     ("jabberwock",
      "ENSEMBL",
      ["name:24", "name:8", "name:28", "nope"],
      {"name:24": ["ENSX26"],
       "name:8": [],
       "name:28": [],
       "nope": []
       }
      )
    ]
)
def test_identifiers_from_symbols_mapping(
        mapper_fixture,
        species,
        authority,
        symbol_list,
        expected_mapping):

    actual = mapping_functions.identifiers_from_symbols_mapping(
        db_path=mapper_fixture.db_path,
        gene_symbol_list=symbol_list,
        species_name=species,
        authority_name=authority
    )
    assert actual['mapping'] == expected_mapping


@pytest.mark.parametrize(
    "species, authority, symbol_list, expected_gene_list",
    [
     ("human",
      "NCBI",
      ["symbol:0", "symbol:6", "symbol:5", "symbol:7", "nope"],
      ["NCBIGene:0", "NCBIGene:6", "test:UNMAPPABLE_NO_MATCH_0",
       "test:UNMAPPABLE_MANY_MATCHES_0",
       "test:UNMAPPABLE_NO_MATCH_1"]
      ),
     ("9606",
      "NCBI",
      ["symbol:0", "symbol:6", "symbol:5", "symbol:7", "nope"],
      ["NCBIGene:0", "NCBIGene:6", "test:UNMAPPABLE_NO_MATCH_0",
       "test:UNMAPPABLE_MANY_MATCHES_0",
       "test:UNMAPPABLE_NO_MATCH_1"]
      ),
     (9606,
      "NCBI",
      ["symbol:0", "symbol:6", "symbol:5", "symbol:7", "nope"],
      ["NCBIGene:0", "NCBIGene:6", "test:UNMAPPABLE_NO_MATCH_0",
       "test:UNMAPPABLE_MANY_MATCHES_0",
       "test:UNMAPPABLE_NO_MATCH_1"]
      ),
     ("jabberwock",
      "ENSEMBL",
      ["symbol:24", "symbol:8", "symbol:28", "nope"],
      ["ENSX26", "test:UNMAPPABLE_NO_MATCH_0", "ENSX30",
       "test:UNMAPPABLE_NO_MATCH_1"]
      ),
     ("jabberwock",
      "ENSEMBL",
      ["name:24", "name:8", "name:28", "nope"],
      ["ENSX26",
       "test:UNMAPPABLE_NO_MATCH_0",
       "test:UNMAPPABLE_NO_MATCH_1",
       "test:UNMAPPABLE_NO_MATCH_2"]
      )
    ]
)
def test_identifiers_from_symbols(
        mapper_fixture,
        species,
        authority,
        symbol_list,
        expected_gene_list):

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        species_obj = query_utils.get_species(
            cursor=cursor,
            species=species
        )

    # now test the method that returns a list of genes
    # with placeholder names applied to unmappable genes

    actual = mapping_functions.identifiers_from_symbols(
        db_path=mapper_fixture.db_path,
        gene_symbol_list=symbol_list,
        species=species_obj,
        authority_name=authority,
        assign_placeholders=True,
        placeholder_prefix="test"
    )
    assert actual['gene_list'] == expected_gene_list


def test_identifiers_from_symbols_mapping_error(
        mapper_fixture):

    # case where there are no citations linking a species
    # to an authority
    msg = "There are 0 citations associated with authority"
    with pytest.raises(ValueError, match=msg):
        mapping_functions.identifiers_from_symbols_mapping(
            db_path=mapper_fixture.db_path,
            gene_symbol_list=["a", "b", "c"],
            species_name="human",
            authority_name="ENSEMBL"
        )

    # case where there is no such species
    msg = "no species match for flotsam"
    with pytest.raises(ValueError, match=msg):
        mapping_functions.identifiers_from_symbols_mapping(
            db_path=mapper_fixture.db_path,
            gene_symbol_list=["a", "b", "c"],
            species_name="flotsam",
            authority_name="NCBI"
        )


def test_identifiers_from_symbols_error(
        mapper_fixture):

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        species = query_utils.get_species(
            cursor=cursor,
            species='human'
        )

    # case where there are no citations linking a species
    # to an authority
    msg = "There are 0 citations associated with authority"
    with pytest.raises(ValueError, match=msg):
        mapping_functions.identifiers_from_symbols(
            db_path=mapper_fixture.db_path,
            gene_symbol_list=["a", "b", "c"],
            species=species,
            authority_name="ENSEMBL"
        )


def test_get_orthologs_mapping_from_ncbi(
        mapper_fixture):

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species='human'
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species='jabberwock'
        )

    gene_idx_list = [0, 1, 2, 4, 7, 6]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapping_functions.ortholog_genes_mapping(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
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

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species=9606
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species='jabberwock'
        )

    actual = mapping_functions.ortholog_genes_mapping(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
        gene_list=gene_list,
        citation_name='NCBI'
    )
    assert actual['mapping'] == expected

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species=9606
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species=999
        )

    actual = mapping_functions.ortholog_genes_mapping(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
        gene_list=gene_list,
        citation_name='NCBI'
    )
    assert actual['mapping'] == expected

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species='mouse'
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species='jabberwock'
        )

    gene_idx_list = [20, 21, 22, 23, 24, 27]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapping_functions.ortholog_genes_mapping(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
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

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species='mouse'
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species='human'
        )

    gene_idx_list = [20, 21, 22, 23, 24, 27]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapping_functions.ortholog_genes_mapping(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
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


@pytest.mark.parametrize(
    "assign_placeholders, placeholder_prefix",
    [(True, None),
     (True, "silly"),
     (False, None)
     ]
)
def test_get_orthologs_from_ncbi(
        mapper_fixture,
        assign_placeholders,
        placeholder_prefix):

    if placeholder_prefix is None:
        prefix = "UNMAPPABLE"
    else:
        prefix = f"{placeholder_prefix}:UNMAPPABLE"

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species='human'
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species='jabberwock'
        )

    gene_idx_list = [0, 1, 2, 4, 7, 6]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapping_functions.ortholog_genes(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
        gene_list=gene_list,
        citation_name='NCBI',
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )

    if assign_placeholders:
        expected = [
            "NCBIGene:14",
            "NCBIGene:13",
            f"{prefix}_NO_MATCH_0",
            "NCBIGene:12",
            "NCBIGene:15",
            f"{prefix}_NO_MATCH_1"
        ]
    else:
        expected = [
            "NCBIGene:14",
            "NCBIGene:13",
            "NCBIGene:2",
            "NCBIGene:12",
            "NCBIGene:15",
            "NCBIGene:6"
        ]

    assert actual['gene_list'] == expected

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species='mouse'
        )
        dst_speices__obj = query_utils.get_species(
            cursor=cursor,
            species='jabberwock'
        )

    gene_idx_list = [20, 21, 22, 23, 24, 27]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapping_functions.ortholog_genes(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
        gene_list=gene_list,
        citation_name='NCBI',
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )

    if assign_placeholders:
        expected = [
            f"{prefix}_NO_MATCH_0",
            "NCBIGene:14",
            f"{prefix}_NO_MATCH_1",
            "NCBIGene:12",
            f"{prefix}_NO_MATCH_2",
            "NCBIGene:13"
        ]
    else:
        expected = [
            "NCBIGene:20",
            "NCBIGene:14",
            "NCBIGene:22",
            "NCBIGene:12",
            "NCBIGene:24",
            "NCBIGene:13"
        ]

    assert actual['gene_list'] == expected

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species='mouse'
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species='human'
        )

    gene_idx_list = [20, 21, 22, 23, 24, 27]
    gene_list = [f'NCBIGene:{ii}' for ii in gene_idx_list]
    actual = mapping_functions.ortholog_genes(
        db_path=mapper_fixture.db_path,
        authority='NCBI',
        src_species=src_species_obj,
        dst_species=dst_species_obj,
        gene_list=gene_list,
        citation_name='NCBI',
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )

    if assign_placeholders:
        expected = [
            f"{prefix}_NO_MATCH_0",
            "NCBIGene:0",
            f"{prefix}_NO_MATCH_1",
            "NCBIGene:4",
            "NCBIGene:6",
            "NCBIGene:1"
        ]
    else:
        expected = [
            "NCBIGene:20",
            "NCBIGene:0",
            "NCBIGene:22",
            "NCBIGene:4",
            "NCBIGene:6",
            "NCBIGene:1"
        ]

    assert actual['gene_list'] == expected


@pytest.mark.parametrize("species_as_int", [True, False])
def test_get_equivalent_genes_mapping_from_ncbi(
        mapper_fixture,
        species_as_int):

    if species_as_int:
        species = "9606"
    else:
        species = "human"

    gene_idx_list = [1, 2, 3, 6, 14, 10]
    gene_list = [f'ENSX{ii}' for ii in gene_idx_list]
    actual = mapping_functions.equivalent_genes_mapping(
        db_path=mapper_fixture.db_path,
        input_authority='ENSEMBL',
        output_authority='NCBI',
        gene_list=gene_list,
        species_name=species,
        citation_name='NCBI'
    )

    expected = {
        'ENSX1': [],
        'ENSX2': ['NCBIGene:1'],
        'ENSX3': [],
        'ENSX6': ['NCBIGene:3'],
        'ENSX14': ['NCBIGene:5', 'NCBIGene:7'],
        'ENSX10': ['NCBIGene:5']
    }

    assert set(expected.keys()) == set(actual['mapping'].keys())
    for k in expected:
        assert set(expected[k]) == set(actual['mapping'][k])


@pytest.mark.parametrize(
    "assign_placeholders, placeholder_prefix, species_as_int",
    [(True, "test", False),
     (True, "test", True),
     (True, "silly", False),
     (False, None, False)]
)
def test_get_equivalent_genes_from_ncbi(
        mapper_fixture,
        assign_placeholders,
        placeholder_prefix,
        species_as_int):

    if species_as_int:
        species = 9606
    else:
        species = "human"

    gene_idx_list = [1, 2, 3, 6, 14, 10]
    gene_list = [f'ENSX{ii}' for ii in gene_idx_list]
    actual = mapping_functions.equivalent_genes(
        db_path=mapper_fixture.db_path,
        input_authority='ENSEMBL',
        output_authority='NCBI',
        gene_list=gene_list,
        species_name=species,
        citation_name='NCBI',
        assign_placeholders=assign_placeholders,
        placeholder_prefix=placeholder_prefix
    )

    if not assign_placeholders:
        expected = [
            "ENSX1",
            "NCBIGene:1",
            "ENSX3",
            "NCBIGene:3",
            "ENSX14",
            "NCBIGene:5"
        ]
    else:
        if placeholder_prefix is not None:
            prefix = f"{placeholder_prefix}:UNMAPPABLE"
        else:
            prefix = "UNMAPPABLE"
        expected = [
            f"{prefix}_NO_MATCH_0",
            "NCBIGene:1",
            f"{prefix}_NO_MATCH_1",
            "NCBIGene:3",
            f"{prefix}_MANY_MATCHES_0",
            "NCBIGene:5"
        ]

    assert actual['gene_list'] == expected


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
def test_alternative_orthologs_mapping(
        mapper_fixture,
        src_species,
        dst_species,
        citation,
        gene_list,
        expected_mapping):

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species=src_species
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species=dst_species
        )

    actual = mapping_functions.ortholog_genes_mapping(
        db_path=mapper_fixture.db_path,
        authority="NCBI",
        src_species=src_species_obj,
        dst_species=dst_species_obj,
        citation_name=citation,
        gene_list=gene_list
    )
    assert actual['mapping'] == expected_mapping


@pytest.mark.parametrize(
    "src_species, dst_species, citation, gene_list, expected_gene_list",
    [("human",
      "jabberwock",
      "alternative_orthologs",
      ["NCBIGene:0", "NCBIGene:2", "NCBIGene:7", "NCBIGene:4", "NCBIGene:6"],
      ["NCBIGene:11", "NCBIGene:13", "NCBIGene:12",
       "UNMAPPABLE_NO_MATCH_0", "UNMAPPABLE_NO_MATCH_1"]
      ),
     ("mouse",
      "human",
      "alternative_orthologs",
      ["NCBIGene:21", "NCBIGene:22", "NCBIGene:23",
       "NCBIGene:24", "NCBIGene:25", "NCBIGene:26",
       "NCBIGene:27"],
      ["UNMAPPABLE_NO_MATCH_0", "UNMAPPABLE_NO_MATCH_1",
       "UNMAPPABLE_NO_MATCH_2", "UNMAPPABLE_NO_MATCH_3",
       "NCBIGene:0", "NCBIGene:2", "UNMAPPABLE_NO_MATCH_4"]
      ),
     ("mouse",
      "human",
      "NCBI",
      ["NCBIGene:21", "NCBIGene:22", "NCBIGene:23",
       "NCBIGene:24", "NCBIGene:25", "NCBIGene:26",
       "NCBIGene:27"],
      ["NCBIGene:0", "UNMAPPABLE_NO_MATCH_0", "NCBIGene:4",
       "NCBIGene:6", "UNMAPPABLE_NO_MATCH_1", "UNMAPPABLE_NO_MATCH_2",
       "NCBIGene:1"]
      ),
     ("mouse",
      "jabberwock",
      "NCBI",
      ["NCBIGene:21", "NCBIGene:22", "NCBIGene:23",
       "NCBIGene:24", "NCBIGene:25", "NCBIGene:26",
       "NCBIGene:27"],
      ["NCBIGene:14", "UNMAPPABLE_NO_MATCH_0", "NCBIGene:12",
       "UNMAPPABLE_NO_MATCH_1", "UNMAPPABLE_NO_MATCH_2",
       "UNMAPPABLE_NO_MATCH_3", "NCBIGene:13"]
      ),
     ]
)
def test_alternative_orthologs(
        mapper_fixture,
        src_species,
        dst_species,
        citation,
        gene_list,
        expected_gene_list):

    with sqlite3.connect(mapper_fixture.db_path) as conn:
        cursor = conn.cursor()
        src_species_obj = query_utils.get_species(
            cursor=cursor,
            species=src_species
        )
        dst_species_obj = query_utils.get_species(
            cursor=cursor,
            species=dst_species
        )

    actual = mapping_functions.ortholog_genes(
        db_path=mapper_fixture.db_path,
        authority="NCBI",
        src_species=src_species_obj,
        dst_species=dst_species_obj,
        citation_name=citation,
        gene_list=gene_list,
        assign_placeholders=True,
        placeholder_prefix=None
    )
    assert actual['gene_list'] == expected_gene_list
