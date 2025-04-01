"""
Define function to just ingest orthologs
"""

import json
import numpy as np
import sqlite3
import time

import mmc_gene_mapper.create_db.utils as db_utils


def insert_orthologs(
        conn,
        gene0_list,
        gene1_list,
        citation_name,
        citation_metadata_dict,
        clobber=False):

    if len(gene0_list) != len(gene1_list):
        raise ValueError(
            f"length of gene lists does not match"
        )

    out_citation = db_utils.insert_unique_citation(
        conn=conn,
        name=citation_name,
        metadata_dict=citation_metadata_dict,
        clobber=clobber
    )

    pair_list = [
        (int(g0), int(g1))
        for g0, g1 in zip(gene0_list, gene1_list)
        if g0 != g1
    ]

    t0 = time.time()
    print(f'=======INGESTING {len(pair_list)} ORTHOLOG PAIRS=======')
    src_citation = db_utils.get_citation(
        conn=conn,
        name='NCBI'
    )

    cursor = conn.cursor()

    # map genes to species
    gene_to_species = dict()
    getter = cursor.execute(
        """
        SELECT
            NCBI_id,
            species_taxon
        FROM NCBI_genes
        WHERE
            citation=?
        """,
        (src_citation['idx'],)
    )

    while True:
        chunk = getter.fetchmany(1000000)
        if len(chunk) == 0:
            break
        for row in chunk:
            if row[0] in gene_to_species:
                if gene_to_species[row[0]] != row[1]:
                    raise ValueError(
                        f"Multiple species for gene {row[0]}"
                    )
            else:
                gene_to_species[row[0]] = row[1]

    print("    GOT GENE TO SPECIES MAP")

    # insert data into ortholog table
    cursor.execute(
        """
        DROP INDEX IF EXISTS ncbi_ortholog_idx
        """
    )

    for ii in range(2):
        if ii == 0:
            i0 = 0
            i1 = 1
        else:
            i0 = 1
            i1 = 0

        values = [
            (pair[i0],
             gene_to_species[pair[i0]],
             pair[i1],
             gene_to_species[pair[i1]],
             out_citation)
            for pair in pair_list
        ]

        cursor.executemany(
            """
            INSERT INTO NCBI_orthologs (
                species0,
                gene0,
                species1,
                gene1,
                citation
            ) VALUES (?, ?, ?, ?, ?)
            """,
            values
        )

    db_utils.create_indexes(conn)
    dur = (time.time()-t0)/60.0
    print(f"=======ORTHOLOG INGESTION TOOK {dur:.2e} minutes=======")
