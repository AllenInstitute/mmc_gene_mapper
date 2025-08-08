"""
Test the utility functions for detecting species from gene symbols
and identifiers
"""

import pytest

import sqlite3

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.mapper.species_detection as species_detection


@pytest.fixture(scope="module")
def symbol_to_species_db_fixture(
        tmp_dir_fixture):
    """
    create a sqlite db file with a gene table and an NCBI_species
    table. This table will be used for testing the function
    that maps gene_symbols to species
    """
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix="symbol_to_species_",
        suffix=".db"
    )

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            "CREATE TABLE NCBI_species (name TEXT, id INT)"
        )
        cursor.execute(
            """
            INSERT INTO NCBI_species (name, id)
            VALUES
                ('mouse', 10090),
                ('human', 9606),
                ('fish', 999),
                ('jabberwock', 9999)
            """
        )

        cursor.execute(
            """
            CREATE TABLE gene
                (symbol TEXT,
                 species_taxon INT,
                 authority INT)
            """
        )
        cursor.execute(
            """
            INSERT INTO gene (symbol, species_taxon, authority)
            VALUES
                ('aaa', 10090, 0),
                ('aaa', 9606, 0),
                ('aaa', 9606, 1),
                ('bbb', 10090, 2),
                ('bbb', 999, 2),
                ('ccc', 999, 2),
                ('ddd', 9606, 3),
                ('eee', 10090, 0),
                ('eee', 9606, 0),
                ('fff', 10090, 4),
                ('fff', 9606, 4)
            """
        )
    return db_path


@pytest.mark.parametrize(
    "gene_list, expected, err_msg",
    [(('aaa', 'bbb', 'eee'),
      {'species': 'mouse', 'species_taxon': 10090},
      None),
     (('aaa', 'bbb', 'eee', 'xxx', 'yyy', 'zzz'),
      {'species': 'mouse', 'species_taxon': 10090},
      None),
     (('aaa', 'bbb', 'eee', 'ddd'),
      None,
      "The gene symbols you gave are consistent with more "
      "than one species"),
     (('xxx', 'yyy', 'zzz'),
      None,
      None)
     ]
)
def test_species_from_symbols(
        symbol_to_species_db_fixture,
        gene_list,
        expected,
        err_msg):

    if err_msg is None:
        actual = species_detection._detect_species_from_symbols(
            db_path=symbol_to_species_db_fixture,
            gene_list=gene_list,
            chunk_size=2
        )
        if expected is not None:
            assert actual == expected
        else:
            assert expected is None
    else:
        with pytest.raises(species_detection.InconsistentSpeciesError,
                           match=err_msg):
            species_detection._detect_species_from_symbols(
                db_path=symbol_to_species_db_fixture,
                gene_list=gene_list,
                chunk_size=2
            )


@pytest.mark.parametrize(
    "gene_list, guess_taxon, expected",
    [(("aaa", "fff", "eee"),
      9606,
      {"species": "human", "species_taxon": 9606}),
     (("aaa", "fff", "eee"),
      10090,
      {"species": "mouse", "species_taxon": 10090}),
     ]
)
def test_species_from_symbols_with_guess(
        symbol_to_species_db_fixture,
        gene_list,
        expected,
        guess_taxon):
    """
    Test that guess_taxon can break a species degeneracy
    in species.
    """
    msg = (
        "Input gene symbols are consistent with several "
        f"species taxons. Taxon {guess_taxon} is one of them. "
        "Will assume the genes are aligned to that taxon."
    )
    with pytest.warns(UserWarning, match=msg):
        actual = species_detection._detect_species_from_symbols(
            db_path=symbol_to_species_db_fixture,
            gene_list=gene_list,
            chunk_size=2,
            guess_taxon=guess_taxon)

    assert actual == expected
