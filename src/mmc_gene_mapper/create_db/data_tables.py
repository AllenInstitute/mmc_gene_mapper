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
            identifier STRING
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
    db_utils.create_index(
        cursor=cursor,
        idx_name="gene_idx",
        table_name="gene",
        column_tuple=(
            "authority",
            "species_taxon",
            "symbol"
        )
    )

    db_utils.create_index(
        cursor=cursor,
        idx_name="gene_equivalence_idx",
        table_name="gene_equivalence",
        column_tuple=(
            "species_taxon",
            "authority0",
            "authority1",
            "gene0"
        )
    )

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
