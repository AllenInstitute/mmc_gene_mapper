"""
Tests to make sure expected errors are thrown when creating
a MMCGeneMapper from a bad file
"""
import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.mapper.mapper as mapper_module


def test_mapper_from_non_file():
    pth = 'not/a/thing'
    msg = "db_path not/a/thing is not a file"
    with pytest.raises(mapper_module.MalformedMapperDBError, match=msg):
        mapper_module.MMCGeneMapper(db_path=pth)


def test_mapper_from_invalid_file(tmp_dir_fixture):
    db_path = file_utils.mkstemp_clean(
       dir=tmp_dir_fixture,
       suffix='.csv'
    )
    with open(db_path, 'w') as dst:
        dst.write('hello')
    with pytest.raises(mapper_module.MalformedMapperDBError):
        mapper_module.MMCGeneMapper(db_path)


def test_mapper_from_invalid_db(tmp_dir_fixture):
    db_path = file_utils.mkstemp_clean(
       dir=tmp_dir_fixture,
       suffix='.db'
    )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute("CREATE TABLE dummy (a INT, b INT)")
        cursor.execute(
            """
            INSERT INTO dummy (a, b) VALUES (4, 5)
            """
        )
    with pytest.raises(mapper_module.MalformedMapperDBError):
        mapper_module.MMCGeneMapper(db_path)


def test_mapper_from_db_with_wrong_validity(tmp_dir_fixture):
    db_path = file_utils.mkstemp_clean(
       dir=tmp_dir_fixture,
       suffix='.db'
    )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "CREATE TABLE mmc_gene_mapper_metadata (validity STR)"
        )
        cursor.execute(
            """
            INSERT INTO mmc_gene_mapper_metadata (validity)
            VALUES ('hello')
            """
        )
    with pytest.raises(mapper_module.MalformedMapperDBError):
        mapper_module.MMCGeneMapper(db_path)
