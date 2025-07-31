"""
This is just a smoke test for the db creation function.

Because this test is expensive and requires downloading data from
ENSEMBL and NCBI, we do not want to run it everytime we run py.test tests/
"""
import pytest

import sqlite3
import tempfile

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.cli.create_db_file as create_db_file


@pytest.fixture(scope='session')
def tmp_dir_fixture(
        tmp_path_factory):
    result = tmp_path_factory.mktemp('mmc_gene_mapper_db_creation_')
    yield result
    print(f'cleaning up {result}')
    file_utils.clean_up(result)


def test_db_creation(tmp_dir_fixture):
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='test_db_',
        suffix='.db'
    )
    local_dir = tempfile.mkdtemp(dir=tmp_dir_fixture)

    create_db_file.create_db_file(
        db_path=db_path,
        local_dir=local_dir,
        ensembl_version=114,
        clobber=True,
        suppress_download_stdout=False,
        n_ensembl=3
    )

    query = """
        SELECT
            gene.identifier,
            gene.symbol
        FROM gene
        WHERE
            gene.authority = ?
        AND
            gene.symbol IS NOT NULL
    """

    # test that the database contains ENSEMBL and NCBI genes
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        for flavor in ("ENSEMBL", "NCBI"):
            print(f"testing {flavor}")
            authority_idx = conn.execute(
                """
                SELECT authority.id
                FROM authority
                WHERE authority.name = ?
                """,
                (flavor,)
            ).fetchall()[0][0]

            result = cursor.execute(query, (authority_idx,))
            chunk = result.fetchmany(10)
            assert len(chunk) > 0
