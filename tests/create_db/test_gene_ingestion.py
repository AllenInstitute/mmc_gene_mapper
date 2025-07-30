import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.data_tables as data_tables
import mmc_gene_mapper.create_db.ncbi_ingestion as ncbi_ingestion


@pytest.fixture(scope='module')
def gene_info_fixture(tmp_dir_fixture):
    """
    Create a gene_info file containing gene symbols
    that pandas may munge from strings into floats.

    Return path to that file
    """
    dst_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='gene_info_',
        suffix='.tsv'
    )
    with open(dst_path, 'w') as dst:
        dst.write('#tax_id\tGeneID\tSymbol\n')
        dst.write('10090\t8\thijkl\n')
        dst.write('10090\t9\tuvwx\n')
        dst.write('10090\t10\tabcd\n')
        dst.write('10090\t11\t7E04\n')
        dst.write('10090\t12\t23.1\n')
    return dst_path


def test_symbols_as_str(
        gene_info_fixture,
        tmp_dir_fixture):
    """
    Test to make sure that gene symbols get ingested from
    gene_info as strings
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='symbols_as_str_',
        suffix='.db'
    )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        data_tables.create_gene_table(cursor)
        ncbi_ingestion.ingest_gene_info(
            conn=conn,
            data_path=gene_info_fixture,
            authority_idx=77,
            citation_idx=88
        )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        results = cursor.execute(
            """
            SELECT
                species_taxon,
                id,
                symbol,
                identifier,
                authority,
                citation
            FROM
                gene
            ORDER BY id
            """
        ).fetchall()

    expected = [
        (10090, 8, 'hijkl', 'NCBIGene:8', 77, 88),
        (10090, 9, 'uvwx', 'NCBIGene:9', 77, 88),
        (10090, 10, 'abcd', 'NCBIGene:10', 77, 88),
        (10090, 11, '7E04', 'NCBIGene:11', 77, 88),
        (10090, 12, '23.1', 'NCBIGene:12', 77, 88)
    ]

    assert results == expected

    for row in results:
        assert isinstance(row[2], str)
        assert isinstance(row[0], int)
        assert isinstance(row[1], int)
