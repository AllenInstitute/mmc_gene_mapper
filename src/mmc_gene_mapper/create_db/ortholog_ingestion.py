"""
Define function to just ingest orthologs
"""

import json
import numpy as np
import pandas as pd
import pathlib
import sqlite3
import time

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.utils as db_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.query_db.query as query_utils


def ingest_hmba_orthologs(
        db_path,
        hmba_file_path,
        citation_name,
        clobber=False,
        chunk_size=1000):

    t0 = time.time()
    db_path = pathlib.Path(db_path)
    if not db_path.is_file():
        raise RuntimeError(
            f"db_path {db_path} is not a file"
        )

    hmba_file_path = pathlib.Path(
        hmba_file_path
    )
    if not hmba_file_path.is_file():
        raise RuntimeError(
            f"{hmba_file_path} is not a file"
        )

    metadata = {
        "file": str(hmba_file_path),
        "hash": file_utils.hash_from_path(hmba_file_path)
    }

    df = pd.read_csv(hmba_file_path)
    gene0_list = df.ncbi_id.values
    gene1_list = df.ortholog_id.values
    tmp_idx_name = "tmp_gene_to_species_idx"
    with sqlite3.connect(db_path) as conn:

        db_utils.create_index(
            cursor=conn.cursor(),
            idx_name=tmp_idx_name,
            table_name="gene",
            column_tuple=("authority", "id")
        )

        ingest_orthologs(
            conn=conn,
            gene0_list=gene0_list,
            gene1_list=gene1_list,
            citation_name=citation_name,
            citation_metadata_dict=metadata,
            clobber=clobber,
            chunk_size=chunk_size
        )

        db_utils.delete_index(
            cursor=conn.cursor(),
            idx_name=tmp_idx_name
        )
    dur = (time.time()-t0)/60.0
    print(f"=======ORTHOLOG INGESTION TOOK {dur:.2e} minutes=======")



def ingest_orthologs(
        conn,
        gene0_list,
        gene1_list,
        citation_name,
        citation_metadata_dict,
        clobber=False,
        chunk_size=1000):

    if len(gene0_list) != len(gene1_list):
        raise ValueError(
            f"length of gene lists does not match"
        )

    out_citation = metadata_utils.insert_unique_citation(
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

    print(f'=======INGESTING {len(pair_list)} ORTHOLOG PAIRS=======')
    src_citation = metadata_utils.get_citation(
        conn=conn,
        name='NCBI'
    )

    src_authority = metadata_utils.get_authority(
        conn=conn,
        name='NCBI'
    )

    cursor = conn.cursor()

    gene_to_species_taxon = dict()
    n_pairs = len(pair_list)
    t0 = time.time()
    n_actual = 0
    for i0 in range(0, n_pairs, chunk_size):
        dur = (time.time()-t0)/60.0
        if i0 > 0:
            per = dur/i0
            pred = per*n_pairs
        else:
            per = 0
            pred = 0
        pair_chunk = pair_list[i0: i0+chunk_size]
        gene_set = set([pair[0] for pair in pair_chunk])
        gene_set = gene_set.union(set([pair[1] for pair in pair_chunk]))
        gene_set = sorted(gene_set - set(gene_to_species_taxon.keys()))

        if len(gene_set) > 0:
            # get gene-to-species map
            query = """
                SELECT
                    id,
                    species_taxon
                FROM gene
                WHERE
                    authority=?
                AND
                    id IN (
            """
            query += ",".join(["?"]*len(gene_set))
            query += ")"
            raw = cursor.execute(
                query,
                (src_authority["idx"],
                 *gene_set)
            ).fetchall()
            for row in raw:
                if row[0] in gene_to_species_taxon:
                    if gene_to_species_taxon[row[0]] != row[1]:
                        raise RuntimeError(
                            "Conflicting species taxon for "
                            f"gene_id {row[0]} authority {authority}; "
                            "unclear how to proceed."
                        )
                gene_to_species_taxon[row[0]] = row[1]


        for order in range(2):
            if order == 0:
                order0 = 0
                order1 = 1
            else:
                order0 = 1
                order1 = 0

            values = [
                (src_authority['idx'],
                 pair[order0],
                 gene_to_species_taxon[pair[order0]],
                 pair[order1],
                 gene_to_species_taxon[pair[order1]],
                 out_citation)
                for pair in pair_chunk
                if pair[order0] in gene_to_species_taxon
                and pair[order1] in gene_to_species_taxon
            ]
            n_actual += len(values)

            cursor.executemany(
                """
                INSERT INTO gene_ortholog (
                    authority,
                    gene0,
                    species0,
                    gene1,
                    species1,
                    citation
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                values
            )


    print(f"=======ACTUALLY INGESTED {n_actual//2} PAIRS=======")
