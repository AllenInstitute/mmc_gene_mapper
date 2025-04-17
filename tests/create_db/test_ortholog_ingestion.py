"""
test ortholog ingestion functions
"""
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.create_db.ortholog_ingestion as ortholog_ingestion


def test_ingest_ortholog_from_species_lookup(
        tmp_dir_fixture):

    # even genes are all connect; odd genes are all connected
    gene0_list = [0, 0, 0, 1, 2, 2, 3, 3]
    gene1_list = [2, 4, 6, 5, 4, 8, 1, 7]

    gene_to_species = {
        ii: ii+4 for ii in range(9)
    }

    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='ortholog_test_db_',
        suffix='.db'
    )

    with sqlite3.connect(db_path) as conn:
        data_utils.create_gene_ortholog_table(conn.cursor())
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
