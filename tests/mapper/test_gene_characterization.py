"""
Tests for function that tries to characterize input gene
identifiers as either ENSEMBL, NCBI, or symbol
"""
import pytest

import numpy as np
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.mapper.gene_characterization as gene_characterization


@pytest.fixture(scope='session')
def gene_characterization_db_fixture(
        tmp_dir_fixture):
    """
    Create a database against which to test gene
    characterization. Return the path to that database.
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix="gene_characteriation_test_",
        suffix=".db"
    )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE authority (id INT, name STRING)
            """
        )
        cursor.execute(
            """
            INSERT INTO authority (id, name)
            VALUES
                (0, 'NCBI'),
                (1, 'ENSEMBL')
            """
        )
        cursor.execute(
            """
            CREATE TABLE gene (
                identifier STRING,
                symbol STRING,
                authority int
            )
            """
        )
        cursor.execute(
            """
            INSERT INTO gene (identifier, symbol, authority)
            VALUES
                ('aaa', 'ENSAP', 0),
                ('ccc', 'ENSAP1', 1),
                ('eee', 'ENSAP2', 0),
                ('ggg', 'NCBIG0', 1),
                ('hhh', 'NCBIG1', 1),
                ('ggg', 'iiii', 0),
                ('jjj', 'aaa', 1)
            """
        )

    return db_path


@pytest.mark.parametrize(
    "gene_list, expected, err_flavor, err_msg",
    [
     (('aaa', 'hhh', 'ENSM999', 'weird',
       'NCBIG0', 'ENSAP1', 'NCBIGene:4', 'unknown',
       'NCBIG1', 'NCBIG88888'),
      ['NCBI', 'ENSEMBL', 'ENSEMBL', 'symbol',
       'symbol', 'symbol', 'NCBI', 'symbol',
       'symbol', 'NCBI'],
      None,
      None),
     (('aaa', 'ccc', 'ggg', 'NCBIG0'),
      None,
      gene_characterization.MultipleAuthorityError,
      "genes mapped to more than one authority"),
    ]
)
def test_dummy(
        gene_characterization_db_fixture,
        gene_list,
        expected,
        err_flavor,
        err_msg):

    if err_flavor is None:
        actual = gene_characterization.characterize_gene_identifiers(
            db_path=gene_characterization_db_fixture,
            gene_list=gene_list,
            chunk_size=4)

        np.testing.assert_array_equal(
            actual=actual,
            desired=np.array(expected)
        )
    else:
        with pytest.raises(err_flavor, match=err_msg):
            gene_characterization.characterize_gene_identifiers(
                db_path=gene_characterization_db_fixture,
                gene_list=gene_list,
                chunk_size=4)
