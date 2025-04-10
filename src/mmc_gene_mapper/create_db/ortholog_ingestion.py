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
import mmc_gene_mapper.create_db.ortholog_utils as ortholog_utils
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
    with sqlite3.connect(db_path) as conn:

        ingest_orthologs_and_citation(
            conn=conn,
            gene0_list=gene0_list,
            gene1_list=gene1_list,
            citation_name=citation_name,
            citation_metadata_dict=metadata,
            clobber=clobber,
            chunk_size=chunk_size
        )

    dur = (time.time()-t0)/60.0
    print(f"=======ORTHOLOG INGESTION TOOK {dur:.2e} minutes=======")


def ingest_orthologs_and_citation(
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

    print(f'=======INGESTING {len(pair_list)} ORTHOLOG PAIRS=======')
    src_citation = metadata_utils.get_citation(
        conn=conn,
        name='NCBI'
    )

    src_authority = metadata_utils.get_authority(
        conn=conn,
        name='NCBI'
    )

    tmp_idx_name = "tmp_gene_to_species_idx"
    db_utils.create_index(
        cursor=conn.cursor(),
        idx_name=tmp_idx_name,
        table_name="gene",
        column_tuple=("authority", "id")
    )

    cursor = conn.cursor()
    gene_to_species_taxon = dict()
    t0 = time.time()
    gene_set = sorted(set(gene0_list).union(set(gene1_list)))

    for i0 in range(0, len(gene_set), chunk_size):
        dur = (time.time()-t0)/60.0
        if i0 > 0:
            per = dur/i0
            pred = per*n_pairs
        else:
            per = 0
            pred = 0
        gene_subset = gene_set[i0:i0+chunk_size]
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
        query += ",".join(["?"]*len(gene_subset))
        query += ")"
        raw = cursor.execute(
            query,
            (src_authority["idx"],
             *gene_subset)
        ).fetchall()
        for row in raw:
            if row[0] in gene_to_species_taxon:
                if gene_to_species_taxon[row[0]] != row[1]:
                    raise ValueError(
                        "Conflicting species taxon for "
                        f"gene_id {row[0]} authority {authority}; "
                        "unclear how to proceed."
                    )
            gene_to_species_taxon[row[0]] = row[1]


    db_utils.delete_index(
        cursor=conn.cursor(),
        idx_name=tmp_idx_name
    )

    _ingest_ortholog_from_species_lookup(
        conn=conn,
        gene0_list=gene0_list,
        gene1_list=gene1_list,
        gene_to_species=gene_to_species_taxon,
        citation_idx=out_citation,
        authority_idx=src_authority['idx']
    )


def ingest_ortholog(
        conn,
        gene0_list,
        gene1_list,
        species0_list,
        species1_list,
        citation_idx,
        authority_idx):
    """
    All lists are numerical indexes
    """

    if len(gene0_list) != len(species0_list):
        raise ValueError(
            "length mismatch between gene0_list and species0_list
        )
    if len(gene1_list) != len(speces1_list):
        raise ValueError(
            "length mismatch between gene1_list and species1_list"
        )

    error_msg = ""
    gene_to_species = dict()
    for gene_list, species_list in [(gene0_list, species0_list),
                                    (gene1_list, species1_list)]"
        for g0, s0 in zip(gene_list, species_list):
            if g0 in gene_to_species:
                if gene_to_species[g0] != s0:
                    error_msg += (
                        f"gene {g0} listed as species {s0} and "
                        f"{gene_to_species[g0]}")
            else:
                gene_to_species[g0] = s0

    if len(error_msg) > 0:
        raise ValueError(error_msg)

    _ingest_ortholog_from_species_lookup(
        conn=conn,
        gene0_list=gene0_list,
        gene1_list=gene1_list,
        gene_to_species=gene_to_species,
        citation_idx=citation_idx,
        authority_idx=authority_idx)


def _ingest_ortholog_from_species_lookup(
        conn,
        gene0_list,
        gene1_list,
        gene_to_species,
        citation_idx,
        authority_idx):

    gene_to_ortholog = ortholog_utils.assig_ortholog_group(
        gene0_list=gene0_list,
        gene1_list=gene1_list
    )

    values = [
        (authority_idx,
         citation_idx,
         gene_to_species[g0],
         g0,
         gene_to_ortholog[g0])
        for g0 in gene_to_species
    ]

    cursor = conn.cursor()
    cursor.executemany(
        """
        INSERT INTO gene_ortholog (
            authority,
            citation,
            species,
            gene,
            ortholog_group
        )
        VALUES (
            ?, ?, ?, ?, ?
        )
        """,
        values
    )
