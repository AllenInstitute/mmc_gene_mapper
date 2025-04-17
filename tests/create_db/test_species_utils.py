import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.create_db.species_utils as species_utils


def test_get_gene_to_species_map(tmp_dir_fixture):

    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='species_map_test_',
        suffix='.db'
    )
    with sqlite3.connect(db_path) as conn:
        data_utils.create_gene_table(conn.cursor())
        cursor = conn.cursor()
        cursor.executemany(
            """
            INSERT INTO gene (
                authority,
                citation,
                id,
                species_taxon
            ) VALUES (?, ?, ?, ?)
            """,
            [
             (0, 0, 0, 0),
             (0, 1, 0, 0),
             (1, 2, 0, 1),
             (0, 0, 1, 1),
             (0, 1, 1, 1),
             (1, 2, 1, 3),
             (0, 0, 2, 0),
             (0, 1, 2, 1)
            ]
        )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        species_map = species_utils.get_gene_to_species_map(
            cursor=cursor,
            gene_list=[0, 1, 55, 66],
            authority_idx=0,
            chunk_size=2
        )
        assert species_map == {
            0: 0,
            1: 1
        }

        species_map = species_utils.get_gene_to_species_map(
            cursor=cursor,
            gene_list=[0, 1, 55, 66],
            authority_idx=1,
            chunk_size=2
        )
        assert species_map == {
            0: 1,
            1: 3
        }

        species_map = species_utils.get_gene_to_species_map(
            cursor=cursor,
            gene_list=[55, 66],
            authority_idx=1,
            chunk_size=2
        )
        assert species_map == dict()

        # test conflicting gene assignments
        msg = "Conflicting species taxon for gene_id 2"
        with pytest.raises(ValueError, match=msg):
            species_utils.get_gene_to_species_map(
                cursor=cursor,
                gene_list=[0, 1, 2, 55, 66],
                authority_idx=0,
                chunk_size=2
            )
