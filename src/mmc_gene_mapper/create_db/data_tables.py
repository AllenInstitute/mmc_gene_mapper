"""
Define utilities around data tables
"""
import mmc_gene_mapper.create_db.utils as db_utils


def create_data_tables(conn):
    print('=======CREATING TABLES=======')
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


def create_data_indexes(conn):
    cursor = conn.cursor()
    create_gene_index(cursor)
    create_gene_equivalence_index(cursor)
    create_gene_ortholog_index(cursor)


def create_gene_index(cursor):
    db_utils.create_index(
        cursor=cursor,
        idx_name="gene_symbol_idx",
        table_name="gene",
        column_tuple=(
            "citation",
            "authority",
            "species_taxon",
            "symbol"
        )
    )

    db_utils.create_index(
        cursor=cursor,
        idx_name="gene_id_idx",
        table_name="gene",
        column_tuple=(
            "citation",
            "authority",
            "species_taxon",
            "id"
        )
    )


def delete_gene_index(cursor):
    db_utils.delete_index(
        cursor=cursor,
        idx_name="gene_symbol_idx"
    )
    db_utils.delete_index(
        cursor=cursor,
        idx_name="gene_id_idx"
    )

def create_gene_equivalence_index(cursor):
    db_utils.create_index(
        cursor=cursor,
        idx_name="gene_equivalence_idx",
        table_name="gene_equivalence",
        column_tuple=(
            "citation",
            "species_taxon",
            "authority0",
            "authority1",
            "gene0"
        )
    )


def delete_gene_equivalence_index(cursor):
    db_utils.delete_index(
        cursor=cursor,
        idx_name="gene_equivalence_idx"
    )


def create_gene_ortholog_index(cursor):
    db_utils.create_index(
        cursor=cursor,
        idx_name="gene_ortholog_idx",
        table_name="gene_ortholog",
        column_tuple=(
            "authority",
            "citation",
            "species0",
            "species1",
            "gene0"
        )
    )

def delete_gene_ortholog_index(cursor):
    db_utils.delete_index(
        cursor=cursor,
        idx_name="gene_ortholog_idx"
    )
