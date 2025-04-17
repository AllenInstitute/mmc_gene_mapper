"""
test ortholog ingestion functions
"""
import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.ortholog_ingestion as ortholog_ingestion


@pytest.mark.parametrize("emit_warning", [True, False])
def test_ingest_ortholog_from_species_lookup(
        tmp_dir_fixture,
        emit_warning):

    # even genes are all connect; odd genes are all connected
    gene0_list = [0, 0, 0, 1, 2, 2, 3, 3]
    gene1_list = [2, 4, 6, 5, 4, 8, 1, 7]

    gene_to_species = {
        ii: ii+4 for ii in range(9)
    }

    if emit_warning:
        gene0_list.append(24)
        gene1_list.append(27)

    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='ortholog_test_db_',
        suffix='.db'
    )

    with sqlite3.connect(db_path) as conn:
        data_utils.create_gene_ortholog_table(conn.cursor())
        if emit_warning:
            msg = "The following genes had no species"
            with pytest.warns(
                    ortholog_ingestion.InvalidOrthologGeneWarning,
                    match=msg):
                ortholog_ingestion._ingest_ortholog_from_species_lookup(
                    conn=conn,
                    gene0_list=gene0_list,
                    gene1_list=gene1_list,
                    gene_to_species=gene_to_species,
                    citation_idx=99,
                    authority_idx=999
                )
        else:
            ortholog_ingestion._ingest_ortholog_from_species_lookup(
                conn=conn,
                gene0_list=gene0_list,
                gene1_list=gene1_list,
                gene_to_species=gene_to_species,
                citation_idx=99,
                authority_idx=999
            )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                authority,
                citation,
                species,
                gene,
                ortholog_group
            FROM gene_ortholog
            ORDER BY gene
            """
        ).fetchall()

    expected = [
        (999, 99, 4, 0, 0),
        (999, 99, 5, 1, 1),
        (999, 99, 6, 2, 0),
        (999, 99, 7, 3, 1),
        (999, 99, 8, 4, 0),
        (999, 99, 9, 5, 1),
        (999, 99, 10, 6, 0),
        (999, 99, 11, 7, 1),
        (999, 99, 12, 8, 0)
    ]
    assert actual == expected


@pytest.fixture
def gene_to_species_fixture():
    """
    Return a dict mapping gene idx to species idx
    (for testing purposes)
    """
    result = {
        ii: ii//3
        for ii in range(12)
    }
    return result


@pytest.fixture
def ortholog_groups_fixture():
    """
    Return lists, each of which is a
    group of orthologs
    """
    return [
        [0, 4, 11],
        [1, 5, 7, 10],
        [2, 9],
        [3, 6]
    ]


def test_ingest_ortholog_from_lists(
        ortholog_groups_fixture,
        gene_to_species_fixture,
        tmp_dir_fixture):
    """
    Test function to ingest orthologs as lists
    of genes and species
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix="ortholog_table_test_",
        suffix=".db"
    )
    gene0_list = [4, 4, 5, 5, 7, 3]
    gene1_list = [0, 11, 1, 7, 10, 6]
    species0_list = [
        gene_to_species_fixture[gene]
        for gene in gene0_list
    ]
    species1_list = [
        gene_to_species_fixture[gene]
        for gene in gene1_list
    ]

    with sqlite3.connect(db_path) as conn:
        data_utils.create_gene_ortholog_table(conn.cursor())
        ortholog_ingestion.ingest_ortholog(
            conn=conn,
            gene0_list=gene0_list,
            gene1_list=gene1_list,
            species0_list=species0_list,
            species1_list=species1_list,
            citation_idx=99,
            authority_idx=999
        )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                authority,
                citation,
                species,
                gene,
                ortholog_group
            FROM gene_ortholog
            ORDER BY gene
            """
        ).fetchall()

    expected = [
        (999, 99, 0, 0, 1),
        (999, 99, 0, 1, 0),
        (999, 99, 1, 3, 2),
        (999, 99, 1, 4, 1),
        (999, 99, 1, 5, 0),
        (999, 99, 2, 6, 2),
        (999, 99, 2, 7, 0),
        (999, 99, 3, 10, 0),
        (999, 99, 3, 11, 1)
    ]

    assert actual == expected


def test_ingest_ortholog_from_lists_errors(
        gene_to_species_fixture,
        tmp_dir_fixture):
    """
    Test that errors are raised if the inputs of ingest_ortholog
    are poorly formed
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix="ortholog_table_test_",
        suffix=".db"
    )
    gene0_list = [4, 4, 5, 5, 7, 3]
    gene0_list_too_big = list(range(8))
    gene1_list = [0, 11, 1, 7, 10, 6]
    gene1_list_too_big = list(range(8))
    species0_list = [
        gene_to_species_fixture[gene]
        for gene in gene0_list
    ]
    species1_list = [
        gene_to_species_fixture[gene]
        for gene in gene1_list
    ]
    species0_list_too_big = [
        gene_to_species_fixture[gene]
        for gene in gene0_list_too_big
    ]

    species0_list_wrong = [0]*len(gene0_list)

    with sqlite3.connect(db_path) as conn:
        data_utils.create_gene_ortholog_table(conn.cursor())
        msg = "length mismatch between gene0_list and species0_list"
        with pytest.raises(ValueError, match=msg):
            ortholog_ingestion.ingest_ortholog(
                conn=conn,
                gene0_list=gene0_list_too_big,
                gene1_list=gene1_list,
                species0_list=species0_list,
                species1_list=species1_list,
                citation_idx=99,
                authority_idx=999
            )

        msg = "length mismatch between gene1_list and species1_list"
        with pytest.raises(ValueError, match=msg):
            ortholog_ingestion.ingest_ortholog(
                conn=conn,
                gene0_list=gene0_list,
                gene1_list=gene1_list_too_big,
                species0_list=species0_list,
                species1_list=species1_list,
                citation_idx=99,
                authority_idx=999
            )

        msg = "length mismatch between gene0_list and gene1_list"
        with pytest.raises(ValueError, match=msg):
            ortholog_ingestion.ingest_ortholog(
                conn=conn,
                gene0_list=gene0_list_too_big,
                gene1_list=gene1_list,
                species0_list=species0_list_too_big,
                species1_list=species1_list,
                citation_idx=99,
                authority_idx=999
            )

        msg = "gene 7 listed as species 2 and 0"
        with pytest.raises(ValueError, match=msg):
            ortholog_ingestion.ingest_ortholog(
                conn=conn,
                gene0_list=gene0_list,
                gene1_list=gene1_list,
                species0_list=species0_list_wrong,
                species1_list=species1_list,
                citation_idx=99,
                authority_idx=999
            )


def test_ingest_orthologs_and_citation(
        gene_to_species_fixture,
        tmp_dir_fixture):
    """
    Test function that ingests orthologs, looking up their species
    from the gene table
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='ortholog_from_gene_test_',
        suffix='.db'
    )

    authority_name='FIAT'

    with sqlite3.connect(db_path) as conn:
        data_utils.create_gene_table(conn.cursor())
        data_utils.create_gene_ortholog_table(conn.cursor())
        metadata_utils.create_metadata_tables(conn)

        # add dummy citations and authorities so that
        # the authority_idx and citation_idx of our test
        # data are non-trivial
        for ii in range(5):
            metadata_utils.insert_authority(
                conn=conn,
                name=f'FAKE{ii}',
                strict=True
            )

        for ii in range(2):
            metadata_utils.insert_citation(
                conn=conn,
                name=f"FAKE_CITATION{ii}",
                metadata_dict={"okay": "fine"}
            )

        authority_idx = metadata_utils.insert_authority(
            conn=conn,
            name=authority_name,
            strict=True
        )
        assert authority_idx == 5

        gene_values = [
            (authority_idx,
             gene_id,
             gene_to_species_fixture[gene_id],
             f'symbol{gene_id}',
             f'identifier:{gene_id}',
             0)
            for gene_id in gene_to_species_fixture.keys()
        ]
        cursor = conn.cursor()
        cursor.executemany(
            """
            INSERT INTO gene (
                authority,
                id,
                species_taxon,
                symbol,
                identifier,
                citation
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            gene_values
        )

    gene0_list = [4, 4, 5, 5, 7, 3]
    gene1_list = [0, 11, 1, 7, 10, 6]

    with sqlite3.connect(db_path) as conn:
        ortholog_ingestion.ingest_orthologs_and_citation(
            conn=conn,
            gene0_list=gene0_list,
            gene1_list=gene1_list,
            citation_name='CITE',
            citation_metadata_dict={'some': 'metdata'},
            clobber=False,
            chunk_size=3,
            gene_authority=authority_name
        )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                authority,
                citation,
                species,
                gene,
                ortholog_group
            FROM gene_ortholog
            ORDER BY gene
            """
        ).fetchall()

    expected = [
        (authority_idx, 2, 0, 0, 1),
        (authority_idx, 2, 0, 1, 0),
        (authority_idx, 2, 1, 3, 2),
        (authority_idx, 2, 1, 4, 1),
        (authority_idx, 2, 1, 5, 0),
        (authority_idx, 2, 2, 6, 2),
        (authority_idx, 2, 2, 7, 0),
        (authority_idx, 2, 3, 10, 0),
        (authority_idx, 2, 3, 11, 1)
    ]
    assert actual == expected
