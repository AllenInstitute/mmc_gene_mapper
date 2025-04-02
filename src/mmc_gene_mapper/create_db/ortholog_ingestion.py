"""
Define function to just ingest orthologs
"""

import json
import numpy as np
import sqlite3
import time

import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.query_db.query as query_utils


def insert_orthologs(
        conn,
        gene0_list,
        species0_list,
        gene1_list,
        species1_list,
        citation_name,
        citation_metadata_dict,
        clobber=False):

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
        (int(g0), int(g1), s0, s1)
        for g0, g1, s0, s1 in zip(
                gene0_list,
                gene1_list,
                species0_list,
                species1_list)
        if g0 != g1
    ]

    t0 = time.time()
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

    species_set = sorted(
        set([p[2] for p in pair_list]).union(
            set([p[3] for p in pair_list])
        )
    )
    species_taxon_lookup = {
        s: query_utils._get_species_taxon(
                cursor=cursor,
                species_name=s,
                strict=True)
        for s in species_set
    }

    print("    GOT SPECIES MAP")

    # insert data into ortholog table
    cursor.execute(
        """
        DROP INDEX IF EXISTS gene_ortholog_idx
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
            (src_authority['idx'],
             pair[i0],
             species_taxon_lookup[pair[i0+2]],
             pair[i1],
             species_taxon_lookup[pair[i1+2]],
             out_citation)
            for pair in pair_list
        ]

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
            values
        )

    data_utils.create_data_indexes(conn)
    dur = (time.time()-t0)/60.0
    print(f"=======ORTHOLOG INGESTION TOOK {dur:.2e} minutes=======")
