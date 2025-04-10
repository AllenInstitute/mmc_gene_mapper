"""
tests for create_db/metadata_tables.py
"""
import pytest

import json
import shutil
import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils


def create_data_tables(conn):
    """
    This is an older version of data_utils.create_data_tables.
    The fact that the schema of gene_ortholog has changed does
    not matter for this unit test
    """
    cursor = conn.cursor()

    cursor.execute(
        """
        CREATE TABLE gene (
            authority INTEGER,
            id INTEGER,
            species_taxon INTEGER,
            symbol STRING,
            identifier STRING,
            citation INTEGER
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE gene_equivalence (
            species_taxon INTEGER,
            authority0 INTEGER,
            gene0 INTEGER,
            authority1 INTEGER,
            gene1 INTEGER,
            citation INTEGER
        )
        """
    )

    cursor.execute(
        """
        CREATE TABLE gene_ortholog(
            authority INTEGER,
            species0 INTEGER,
            gene0 INTEGER,
            species1 INTEGER,
            gene1 INTEGER,
            citation INTEGER
        )
        """
    )


@pytest.fixture(scope='module')
def citation_dataset_fixture():
    """
    Data to be ingested into test citation table
    """
    return [
        {"name": "A",
         "metadata": {"A": 1, "B": {"x": "a", "y": 2}},
         "idx": 0},
        {"name": "B",
         "metadata": {"some": "metadata"},
         "idx": 1},
        {"name": "C",
         "metadata": {"other": "metadata"},
         "idx": 2},
        {"name": "C",
         "metadata": {"uh": "oh"},
         "idx": 3}
    ]


@pytest.fixture(scope='module')
def authority_dataset_fixture():
    """
    Data to be ingested into authority table
    """
    return [
        {"name": "a0", "idx": 0},
        {"name": "a1", "idx": 1},
        {"name": "a2", "idx": 2},
        {"name": "a2", "idx": 3}
    ]


@pytest.fixture(scope='module')
def gene_data_fixture():
    """
    List of records for gene table that reference
    simulated citations

    column order = (
        authority
        id
        species_taxon
        symbol
        identifier
        citation
    )
    """
    return [
        (0, 1, 11, 'ENS111', 111, 0),
        (0, 2, 22, 'ENS222', 222, 0),
        (1, 3, 33, 'ENS333', 333, 1)
    ]


@pytest.fixture(scope='module')
def gene_equivalence_data_fixture():
    """
    List of records for gene_equivalence table
    that reference simulated citations

    column order = (
        species_taxon
        authority0
        gene0
        authority1
        gene1
        citation
    )
    """
    return [
        (11, 0, 1, 1, 2, 0),
        (11, 1, 2, 2, 3, 0),
        (22, 2, 4, 0, 5, 1),
        (33, 2, 6, 3, 7, 2)
    ]


@pytest.fixture(scope='module')
def gene_ortholog_data_fixture():
    """
    List of records for gene_ortholog table
    that reference simulated citations

    column_order= (
        authority
        species0
        gene0
        species1
        gene1
        citation
    )
    """
    return [
        (0, 1, 1, 2, 2, 0),
        (1, 3, 3, 4, 4, 1),
        (1, 5, 5, 6, 6, 0),
        (0, 7, 7, 8, 8, 2)
    ]


@pytest.fixture(scope='module')
def metadata_table_fixture(
        citation_dataset_fixture,
        authority_dataset_fixture,
        tmp_dir_fixture):
    """
    Return path to sqlite database with citation
    table and authority table to test against.
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='citation_table_',
        suffix='.db'
    )
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE citation (
                id INTEGER,
                name STRING,
                metadata STRING
            )
            """
        )
        values = [
            (row['idx'], row['name'], json.dumps(row['metadata']))
            for row in citation_dataset_fixture
        ]
        cursor.executemany(
            """
            INSERT INTO citation (
                id,
                name,
                metadata
            ) VALUES (?, ?, ?)
            """,
            values
        )

        cursor.execute(
            """
            CREATE TABLE authority (
                id INTEGER,
                name STRING
            )
            """
        )

        cursor.executemany(
            """
            INSERT INTO authority (
                id,
                name
            ) VALUES (?, ?)
            """,
            [(row["idx"], row["name"])
             for row in authority_dataset_fixture]
        )

    return db_path


@pytest.fixture(scope='function')
def metadata_table_with_data_fixture(
        metadata_table_fixture,
        gene_data_fixture,
        gene_equivalence_data_fixture,
        gene_ortholog_data_fixture,
        tmp_dir_fixture):
    """
    Return path to a database that has a citation
    table and data tables.

    scope is 'function' because we are expected to change it
    in tests
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='citation_with_data_',
        suffix='.db'
    )
    shutil.copy(
        src=metadata_table_fixture,
        dst=db_path
    )

    with sqlite3.connect(db_path) as conn:
        create_data_tables(conn)
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
            gene_data_fixture
        )

        cursor.executemany(
            """
            INSERT INTO gene_equivalence (
                species_taxon,
                authority0,
                gene0,
                authority1,
                gene1,
                citation
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            gene_equivalence_data_fixture
        )

        cursor.executemany(
            """
            INSERT INTO gene_ortholog (
                authority,
                species0,
                gene0,
                species1,
                gene1,
                citation
            ) VALUES (?, ?, ?, ?, ?, ?)
            """,
            gene_ortholog_data_fixture
        )

    return db_path


@pytest.mark.parametrize("strict", [True, False])
def test_get_citation(
        citation_dataset_fixture,
        metadata_table_fixture,
        strict):

    with sqlite3.connect(metadata_table_fixture) as conn:
        a_result = metadata_utils.get_citation(
            conn=conn,
            name="A",
            strict=strict)
        assert a_result == citation_dataset_fixture[0]

        a_result = metadata_utils.get_citation(
            conn=conn,
            name="B",
            strict=strict)
        assert a_result == citation_dataset_fixture[1]

        with pytest.raises(ValueError, match="More than one citation"):
            metadata_utils.get_citation(
                conn=conn,
                name="C",
                strict=strict
            )

        if strict:
            with pytest.raises(ValueError, match="No citation"):
                metadata_utils.get_citation(
                    conn=conn,
                    name="D",
                    strict=strict
                )
        else:
            assert metadata_utils.get_citation(
                conn=conn,
                name="D",
                strict=strict) is None


def test_delete_citation(
        metadata_table_with_data_fixture,
        citation_dataset_fixture,
        gene_data_fixture,
        gene_equivalence_data_fixture,
        gene_ortholog_data_fixture):
    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        metadata_utils.delete_citation(
            conn=conn,
            name="A"
        )

    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        cursor = conn.cursor()
        citations = cursor.execute(
            "SELECT name, metadata, id FROM citation"
        ).fetchall()
        assert len(citations) == 3
        expected = set([
            (row['name'],
             json.dumps(row['metadata']),
             row['idx'])
            for row in citation_dataset_fixture[1:]
        ])
        assert set(citations) == expected

        genes = cursor.execute(
            """
            SELECT
                authority,
                id,
                species_taxon,
                symbol,
                identifier,
                citation
            FROM
                gene
            """
        ).fetchall()
        assert set(genes) == set(gene_data_fixture[2:])

        equiv = cursor.execute(
            """
            SELECT
                species_taxon,
                authority0,
                gene0,
                authority1,
                gene1,
                citation
            FROM
                gene_equivalence
            """
        ).fetchall()
        assert set(equiv) == set(gene_equivalence_data_fixture[2:])

        orthologs = cursor.execute(
            """
            SELECT
                authority,
                species0,
                gene0,
                species1,
                gene1,
                citation
            FROM
                gene_ortholog
            """
        ).fetchall()
        assert set(orthologs) == set([
            gene_ortholog_data_fixture[1],
            gene_ortholog_data_fixture[3]])


def test_insert_citation(
        metadata_table_with_data_fixture,
        citation_dataset_fixture):
    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        with pytest.raises(ValueError, match="already exists"):
            metadata_utils.insert_citation(
                conn=conn,
                name="A",
                metadata_dict={"will not": "work"}
            )

        assert metadata_utils.insert_citation(
            conn=conn,
            name="E",
            metadata_dict={"okay": "fine"}
        ) == 4

        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                name,
                metadata,
                id
            FROM
                citation
            """
        ).fetchall()

        expected = [
            (row['name'],
             json.dumps(row['metadata']),
             row['idx'])
            for row in citation_dataset_fixture
        ]
        expected.append(
            ("E",
             json.dumps({"okay": "fine"}),
             4)
        )
        assert len(expected) == len(actual)
        for ex in expected:
            assert ex in actual


@pytest.mark.parametrize("clobber", [True, False])
def test_insert_unique_citation(
        metadata_table_with_data_fixture,
        citation_dataset_fixture,
        gene_data_fixture,
        gene_equivalence_data_fixture,
        gene_ortholog_data_fixture,
        clobber):
    """
    Test that if we insert a unique citation with name='A'
    and clobber=True, all of the relevant data gets deleted
    and the new 'A' citation exists
    """
    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        if not clobber:
            msg = "already exists; run with clobber"
            with pytest.raises(ValueError, match=msg):
                metadata_utils.insert_unique_citation(
                    conn=conn,
                    name="A",
                    metadata_dict={"alt": "A"},
                    clobber=clobber
                )
        else:
            metadata_utils.insert_unique_citation(
                conn=conn,
                name="A",
                metadata_dict={"alt": "A"},
                clobber=clobber
            )

    if not clobber:
        return

    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        cursor = conn.cursor()
        citations = cursor.execute(
            "SELECT name, metadata, id FROM citation"
        ).fetchall()
        assert len(citations) == 4
        expected = set([
            (row['name'],
             json.dumps(row['metadata']),
             row['idx'])
            for row in citation_dataset_fixture[1:]
        ])
        expected.add(
            ("A",
             json.dumps({"alt": "A"}),
             4)
        )
        assert set(citations) == expected

        genes = cursor.execute(
            """
            SELECT
                authority,
                id,
                species_taxon,
                symbol,
                identifier,
                citation
            FROM
                gene
            """
        ).fetchall()
        assert set(genes) == set(gene_data_fixture[2:])

        equiv = cursor.execute(
            """
            SELECT
                species_taxon,
                authority0,
                gene0,
                authority1,
                gene1,
                citation
            FROM
                gene_equivalence
            """
        ).fetchall()
        assert set(equiv) == set(gene_equivalence_data_fixture[2:])

        orthologs = cursor.execute(
            """
            SELECT
                authority,
                species0,
                gene0,
                species1,
                gene1,
                citation
            FROM
                gene_ortholog
            """
        ).fetchall()
        assert set(orthologs) == set([
            gene_ortholog_data_fixture[1],
            gene_ortholog_data_fixture[3]])


@pytest.mark.parametrize("strict", [True, False])
def test_get_authority(
        metadata_table_fixture,
        strict):

    with sqlite3.connect(metadata_table_fixture) as conn:

        aa = metadata_utils.get_authority(
            conn=conn,
            name="a1",
            strict=strict)

        assert aa == {"name": "a1", "idx": 1}

        with pytest.raises(ValueError, match="More than one authority"):
            metadata_utils.get_authority(
                conn=conn,
                name="a2",
                strict=strict
            )

        if strict:
            with pytest.raises(ValueError, match="No authority"):
                metadata_utils.get_authority(
                    conn=conn,
                    name="b0",
                    strict=strict
                )
        else:
            assert metadata_utils.get_authority(
                conn=conn,
                name="b0",
                strict=strict
            ) is None


def test_delete_authority(
        metadata_table_with_data_fixture,
        authority_dataset_fixture,
        gene_data_fixture,
        gene_equivalence_data_fixture,
        gene_ortholog_data_fixture):

    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        metadata_utils.delete_authority(
            conn=conn,
            name="a0"
        )

    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        cursor = conn.cursor()
        citations = cursor.execute(
            "SELECT name, id FROM authority"
        ).fetchall()
        assert len(citations) == 3
        expected = set([
            (row['name'],
             row['idx'])
            for row in authority_dataset_fixture[1:]
        ])
        assert set(citations) == expected

        genes = cursor.execute(
            """
            SELECT
                authority,
                id,
                species_taxon,
                symbol,
                identifier,
                citation
            FROM
                gene
            """
        ).fetchall()
        assert set(genes) == set(gene_data_fixture[2:])

        equiv = cursor.execute(
            """
            SELECT
                species_taxon,
                authority0,
                gene0,
                authority1,
                gene1,
                citation
            FROM
                gene_equivalence
            """
        ).fetchall()
        assert set(equiv) == set(
            [gene_equivalence_data_fixture[1],
             gene_equivalence_data_fixture[3]])

        orthologs = cursor.execute(
            """
            SELECT
                authority,
                species0,
                gene0,
                species1,
                gene1,
                citation
            FROM
                gene_ortholog
            """
        ).fetchall()
        assert set(orthologs) == set([
            gene_ortholog_data_fixture[1],
            gene_ortholog_data_fixture[2]])


@pytest.mark.parametrize("strict", [True, False])
def test_insert_authority(
        metadata_table_with_data_fixture,
        authority_dataset_fixture,
        strict):
    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        if strict:
            with pytest.raises(ValueError, match="already exists"):
                metadata_utils.insert_authority(
                    conn=conn,
                    name="a1",
                    strict=strict
                )
        else:
            assert metadata_utils.insert_authority(
                conn=conn,
                name="a1",
                strict=strict) == 1

        assert metadata_utils.insert_authority(
            conn=conn,
            name="c0",
            strict=strict
        ) == 4

        cursor = conn.cursor()
        actual = cursor.execute(
            """
            SELECT
                name,
                id
            FROM
                authority
            """
        ).fetchall()

        expected = [
            (row['name'],
             row['idx'])
            for row in authority_dataset_fixture
        ]
        expected.append(
            ("c0",
             4)
        )
        assert len(expected) == len(actual)
        for ex in expected:
            assert ex in actual


@pytest.mark.parametrize("clobber", [True, False])
def test_insert_unique_authority(
        metadata_table_with_data_fixture,
        authority_dataset_fixture,
        gene_data_fixture,
        gene_equivalence_data_fixture,
        gene_ortholog_data_fixture,
        clobber):
    """
    Test that if we insert a unique authority with name='a0'
    and clobber=True, all of the relevant data gets deleted
    and the new 'a0' authority exists
    """
    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        if not clobber:
            msg = "already exists; run with clobber"
            with pytest.raises(ValueError, match=msg):
                metadata_utils.insert_unique_authority(
                    conn=conn,
                    name="a0",
                    clobber=clobber
                )
        else:
            metadata_utils.insert_unique_authority(
                conn=conn,
                name="a0",
                clobber=clobber
            )

    if not clobber:
        return

    with sqlite3.connect(metadata_table_with_data_fixture) as conn:
        cursor = conn.cursor()
        citations = cursor.execute(
            "SELECT name, id FROM authority"
        ).fetchall()
        assert len(citations) == 4
        expected = set([
            (row['name'],
             row['idx'])
            for row in authority_dataset_fixture[1:]
        ])
        expected.add(
            ("a0",
             4)
        )
        assert set(citations) == expected

        genes = cursor.execute(
            """
            SELECT
                authority,
                id,
                species_taxon,
                symbol,
                identifier,
                citation
            FROM
                gene
            """
        ).fetchall()
        assert set(genes) == set(gene_data_fixture[2:])

        equiv = cursor.execute(
            """
            SELECT
                species_taxon,
                authority0,
                gene0,
                authority1,
                gene1,
                citation
            FROM
                gene_equivalence
            """
        ).fetchall()
        assert set(equiv) == set([
            gene_equivalence_data_fixture[1],
            gene_equivalence_data_fixture[3]
        ])

        orthologs = cursor.execute(
            """
            SELECT
                authority,
                species0,
                gene0,
                species1,
                gene1,
                citation
            FROM
                gene_ortholog
            """
        ).fetchall()
        assert set(orthologs) == set([
            gene_ortholog_data_fixture[1],
            gene_ortholog_data_fixture[2]])
