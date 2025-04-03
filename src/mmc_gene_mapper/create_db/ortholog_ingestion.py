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
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.data_tables as data_utils
import mmc_gene_mapper.query_db.query as query_utils


def ingest_hmba_orthologs(
        db_path,
        hmba_file_path,
        citation_name,
        baseline_species='human',
        clobber=False):

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
    species0_list=[
        s.replace('_', ' ')
        for s in df.species.values
    ]

    gene1_list = df.ortholog_id.values
    species1_list=[baseline_species]*len(gene1_list)
    with sqlite3.connect(db_path) as conn:
        ingest_orthologs(
            conn=conn,
            gene0_list=gene0_list,
            species0_list=species0_list,
            gene1_list=gene1_list,
            species1_list=species1_list,
            citation_name=citation_name,
            citation_metadata_dict=metadata,
            clobber=clobber
        )


def ingest_orthologs(
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

    dur = (time.time()-t0)/60.0
    print(f"=======ORTHOLOG INGESTION TOOK {dur:.2e} minutes=======")
