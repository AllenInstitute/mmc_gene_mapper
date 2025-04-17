"""
Define function to just ingest orthologs
"""

import pandas as pd
import pathlib
import sqlite3
import time
import warnings

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.utils as db_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.create_db.ortholog_utils as ortholog_utils
import mmc_gene_mapper.create_db.species_utils as species_utils


def ingest_hmba_orthologs(
        db_path,
        hmba_file_path,
        citation_name,
        primary_id_column='ncbi_id',
        ortholog_id_column='ortholog_id',
        gene_authority='NCBI',
        clobber=False,
        chunk_size=1000):
    """
    Ingest orthologs from a CSV file containing a primary_id
    column and a ortholog_id column

    Parameters
    ----------
    db_path:
        path to the database into which we are ingesting
    hmba_file_path:
        path to the CSV containing the ortholog information
    citation_name:
        the name that will be given to the citation associated
        with this file
    primary_id_column:
        column in the CSV file giving the integer ID of the
        "anchor" genes
    ortholog_id_column:
        column in the CSV file giving the integer ID of the
        genes that are orthologs to the genes in primary_id_column
    gene_authority:
        the human-readable name of the authority in which these gene
        IDs are defined
    clobber:
        if True and a citation with this citation_name already exists,
        overwrite; raise an exception otherwise
    chunk_size:
        the number of genes to ingest at a time (to avoid overwhelming
        machine memory in the case of large files)

    Returns
    -------
    None
        orthologs are populated in gene_ortholog table of the specified
        database file
    """

    hmba_file_path = pathlib.Path(hmba_file_path)
    print(f'=======INGESTING ORTHOLOGS FROM {hmba_file_path.name}=======')

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
    df = df[df[primary_id_column] != df[ortholog_id_column]]
    gene0_list = [int(ii) for ii in df[primary_id_column].values]
    gene1_list = [int(ii) for ii in df[ortholog_id_column].values]
    with sqlite3.connect(db_path) as conn:

        ingest_orthologs_creating_citation(
            conn=conn,
            gene0_list=gene0_list,
            gene1_list=gene1_list,
            citation_name=citation_name,
            citation_metadata_dict=metadata,
            clobber=clobber,
            chunk_size=chunk_size,
            gene_authority=gene_authority
        )

    dur = (time.time()-t0)/60.0
    print(f"=======ORTHOLOG INGESTION TOOK {dur:.2e} minutes=======")


def ingest_ncbi_orthologs(
        conn,
        data_path,
        citation_idx,
        authority_idx):
    """
    Ingest a file that looks like NCBI's gene_orthologs file.

    Parameters
    ----------
    conn:
        sqlite3 connection to the database
    data_path:
        path to the gene_orthologs file being ingested
    citation_idx:
        the integer identifying the citation that is to be associated
        with this set of orthologs
    authority_idx:
        the integer identifying the authority that is to be associated
        with this set of orthologs

    Returns
    -------
    None
        orthologs are ingested into the gene_ortholog table of the database
    """
    t0 = time.time()
    print('=======INGESTING ORTHOLOGS=======')

    data = pd.read_csv(data_path, delimiter='\t')
    data = data[data['relationship'] == 'Ortholog']

    gene0_list = [
        int(ii) for ii in data['GeneID'].values
    ]
    gene1_list = [
        int(ii) for ii in data['Other_GeneID'].values
    ]

    ingest_orthologs_specifying_citation(
        conn=conn,
        gene0_list=gene0_list,
        gene1_list=gene1_list,
        citation_idx=citation_idx,
        authority_idx=authority_idx
    )

    dur = (time.time()-t0)/60.0
    print(f'=======INGESTING gene_orthologs TOOK {dur:.2e} minutes=======')


def ingest_orthologs_creating_citation(
        conn,
        gene0_list,
        gene1_list,
        citation_name,
        citation_metadata_dict,
        clobber=False,
        chunk_size=1000,
        gene_authority='NCBI'):
    """
    Ingest orthologs, creating a citation for this batch of
    orthologs.

    Parameters
    ----------
    conn:
        sqlite3 connection object
    gene0_list:
        list of integers identifying genes
    gene1_list:
        list of integers identifying genes that are orthologs
        of genes in gene0_list (gene1_list[ii] is ortholog of
        gene0_list[ii])
    citation_name:
        name to be assigned to the citation created for these
        orthologs
    citation_metadata_dict:
        dict defining the metadata to be associated with the
        citation created for these orthologs
    clobber:
        if False and a citation with this name already exists,
        fail. If not, overwrite the existing citation (deleting
        all data connected with that citation)
    chunk_size:
        genes to process at once (to avoid overwhelming machine
        memory)
    gene_authority:
        the human-readable name of the authority in which
        these genes are defined.

    """

    if len(gene0_list) != len(gene1_list):
        raise ValueError(
            "length of gene lists does not match"
        )

    out_citation = metadata_utils.insert_unique_citation(
        conn=conn,
        name=citation_name,
        metadata_dict=citation_metadata_dict,
        clobber=clobber
    )

    src_authority = metadata_utils.get_authority(
        conn=conn,
        name=gene_authority,
        strict=True
    )

    ingest_orthologs_specifying_citation(
        conn=conn,
        gene0_list=gene0_list,
        gene1_list=gene1_list,
        citation_idx=out_citation,
        authority_idx=src_authority['idx'],
        chunk_size=chunk_size
    )


def ingest_orthologs_specifying_citation(
        conn,
        gene0_list,
        gene1_list,
        citation_idx,
        authority_idx,
        chunk_size=10000):
    """
    Ingest ortholog data into database

    Prameters
    ---------
    conn:
        the sqlite3 connection to the database
    gene0_list:
        list of integers; the first set of genes
    gene1_list:
        list of integers; genes that are orthologs to the
        corresponding genes in gene0_list
    citation_idx:
        the integer corresponding to the citation justifying
        these ortholog assignments
    authority_idx:
        the integer corresponding to the authority in which
        these orthologs are identified (ENSEMBL or NCBI)
    chunk_size:
        genes to process at once (to avoid overwhelming machine
        memory)

    Returns
    -------
    None
        data is ingested into the gene_ortholog database
    """
    tmp_idx_name = "tmp_gene_to_species_idx"
    db_utils.create_index(
        cursor=conn.cursor(),
        idx_name=tmp_idx_name,
        table_name="gene",
        column_tuple=("authority", "id")
    )

    cursor = conn.cursor()
    gene_set = sorted(set(gene0_list).union(set(gene1_list)))

    gene_to_species_taxon = species_utils.get_gene_to_species_map(
        cursor=cursor,
        gene_list=gene_set,
        authority_idx=authority_idx,
        chunk_size=chunk_size
    )

    db_utils.delete_index(
        cursor=conn.cursor(),
        idx_name=tmp_idx_name
    )

    _ingest_orthologs_from_species_lookup(
        conn=conn,
        gene0_list=gene0_list,
        gene1_list=gene1_list,
        gene_to_species=gene_to_species_taxon,
        citation_idx=citation_idx,
        authority_idx=authority_idx
    )


def _ingest_orthologs_from_species_list(
        conn,
        gene0_list,
        gene1_list,
        species0_list,
        species1_list,
        citation_idx,
        authority_idx):
    """
    Ingest ortholog data into database

    Prameters
    ---------
    conn:
        the sqlite3 connection to the database
    gene0_list:
        list of integers; the first set of genes
    gene1_list:
        list of integers; genes that are orthologs to the
        corresponding genes in gene0_list
    species0_list:
        list of integers indicating the species that the
        genes in gene0_list belong to
    species1_list:
        list of integers indicating the specie that the
        genes in gene1_list belong to
    citation_idx:
        the integer corresponding to the citation justifying
        these ortholog assignments
    authority_idx:
        the integer corresponding to the authority in which
        these orthologs are identified (ENSEMBL or NCBI)

    Returns
    -------
    None
        data is ingested into the gene_ortholog database
    """

    if len(gene0_list) != len(species0_list):
        raise ValueError(
            "length mismatch between gene0_list and species0_list"
        )
    if len(gene1_list) != len(species1_list):
        raise ValueError(
            "length mismatch between gene1_list and species1_list"
        )
    if len(gene1_list) != len(gene0_list):
        raise ValueError(
            "length mismatch between gene0_list and gene1_list"
        )

    error_msg = ""
    gene_to_species = dict()
    for gene_list, species_list in [(gene0_list, species0_list),
                                    (gene1_list, species1_list)]:
        for g0, s0 in zip(gene_list, species_list):
            if g0 in gene_to_species:
                if gene_to_species[g0] != s0:
                    error_msg += (
                        f"gene {g0} listed as species {s0} and "
                        f"{gene_to_species[g0]}")
            else:
                gene_to_species[g0] = s0

    print("    GOT gene_to_species")

    if len(error_msg) > 0:
        raise ValueError(error_msg)

    _ingest_orthologs_from_species_lookup(
        conn=conn,
        gene0_list=gene0_list,
        gene1_list=gene1_list,
        gene_to_species=gene_to_species,
        citation_idx=citation_idx,
        authority_idx=authority_idx)


def _ingest_orthologs_from_species_lookup(
        conn,
        gene0_list,
        gene1_list,
        gene_to_species,
        citation_idx,
        authority_idx):
    """
    Ingest ortholog data into database

    Prameters
    ---------
    conn:
        the sqlite3 connection to the database
    gene0_list:
        list of integers; the first set of genes
    gene1_list:
        list of integers; genes that are orthologs to the
        corresponding genes in gene0_list
    gene_to_species:
        a dict mapping the gene identifier ints in gene0_list
        and gene1_list to species identifier ints
    citation_idx:
        the integer corresponding to the citation justifying
        these ortholog assignments
    authority_idx:
        the integer corresponding to the authority in which
        these orthologs are identified (ENSEMBL or NCBI)

    Returns
    -------
    None
        data is ingested into the gene_ortholog database

    Notes
    -----
    genes that do not occur in gene_to_species will not be ingested
    """

    print(f'gene0 {len(gene0_list)}')
    print(f'gene1 {len(gene1_list)}')

    gene_to_ortholog = ortholog_utils.assign_ortholog_group(
        gene0_list=gene0_list,
        gene1_list=gene1_list
    )

    print(f'gene_to_ortholog {len(gene_to_ortholog)}')
    print(f'gene_to_species {len(gene_to_species)}')

    warning_msg = ""
    for gene in gene_to_ortholog:
        if gene not in gene_to_species:
            warning_msg += f"{gene}\n"
    if len(warning_msg) > 0:
        warning_msg = (
            "The following genes had no species assigned to "
            "them and were not ingested\n"
            f"{warning_msg}"
        )
        warnings.warn(warning_msg, category=InvalidOrthologGeneWarning)

    values = [
        (authority_idx,
         citation_idx,
         gene_to_species[g0],
         g0,
         gene_to_ortholog[g0])
        for g0 in gene_to_species
    ]

    print(f"inserting {len(values)} values with citation {citation_idx}")
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


class InvalidOrthologGeneWarning(UserWarning):
    pass
