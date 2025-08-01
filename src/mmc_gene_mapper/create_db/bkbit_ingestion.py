"""
Define function to ingest genes from bkbit file
"""

import json
import pathlib
import sqlite3

import mmc_gene_mapper.utils.str_utils as str_utils
import mmc_gene_mapper.create_db.utils as db_utils
import mmc_gene_mapper.create_db.metadata_tables as metadata_utils
import mmc_gene_mapper.query_db.query as query_utils


def ingest_bkbit_genes(
        db_path,
        bkbit_path):
    """
    Ingest gene info from a bkbit jsonld file into the gene
    table of the database

    Parameters
    ----------
    db_path:
        path to the database into which we are ingesting
    bkbit_path:
        path to the jsonld file with the data being ingested

    Returns
    -------
    None
    """

    bkbit_path = pathlib.Path(bkbit_path)
    if not bkbit_path.is_file():
        raise RuntimeError(
            f"{bkbit_path} is not a file"
        )
    db_path = pathlib.Path(db_path)
    if not db_utils.check_existence(db_path):
        raise RuntimeError(
            f"{db_path} does not exist"
        )

    metadata_dict, values, citation_name = read_bkbit_data(
        bkbit_path=bkbit_path,
        db_path=db_path
    )

    with sqlite3.connect(db_path) as conn:
        citation_idx = metadata_utils.insert_unique_citation(
            conn=conn,
            name=citation_name,
            metadata_dict=metadata_dict,
            clobber=False
        )
        values = [
            (*r, citation_idx) for r in values
        ]
        cursor = conn.cursor()
        cursor.execute("DROP INDEX IF EXISTS gene_idx")
        cursor.executemany(
            """
            INSERT INTO gene (
                authority,
                id,
                species_taxon,
                symbol,
                identifier,
                citation
            )
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            values
        )


def read_bkbit_data(bkbit_path, db_path):
    """
    Read data from a bkbit file and format for ingestion into
    the database.

    Parameters
    ----------
    bkbit_path:
        path to the bkbit jsonld file
    db_path:
        path to the database into which we are ingesting data
        (for looking up the mapping from species name to
        species_taxon)

    Returns
    -------
    A tuple
        (
            metadata -- dict of metadata for the citation
                corresponding to this file
            raw_gene_annotation_values -- list of tuples; the
                values to be ingested into the gene table
            citation_name -- the human readable name of the
                citation corresponding to this file
        )
    """
    with open(bkbit_path, "rb") as src:
        bkbit_data = json.load(src)

    metadata_dict = dict()
    authority_idx = None
    species_lookup = dict()
    raw_gene_annotation_values = []
    citation_name = None

    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        for record in bkbit_data['@graph']:
            if 'bican:Checksum' in record['category']:
                continue
            elif 'bican:GeneAnnotation' in record['category']:

                if 'symbol' not in record:
                    record['symbol'] = None
                if 'name' not in record:
                    record['name'] = None

                try:
                    gene_idx = str_utils.int_from_identifier(
                        record['source_id']
                    )
                except str_utils.MalformedGeneIdentifierError:
                    continue

                if record['in_taxon_label'] not in species_lookup:
                    species_lookup[record['in_taxon_label']] = (
                        query_utils._get_species_taxon(
                            cursor=cursor,
                            species_name=record['in_taxon_label']
                        )
                    )
                species_idx = species_lookup[record['in_taxon_label']]
                if authority_idx is None:
                    raise RuntimeError(
                        "authority_idx is None"
                    )

                gene = (
                    authority_idx,
                    gene_idx,
                    species_idx,
                    record['symbol'],
                    record['source_id']
                )
                raw_gene_annotation_values.append(gene)
                if record['name'] != record['symbol']:
                    gene = (
                        authority_idx,
                        gene_idx,
                        species_idx,
                        record['name'],
                        record['source_id']
                    )
                    raw_gene_annotation_values.append(gene)
            else:
                category = record['category']
                if len(category) > 1:
                    raise RuntimeError(
                        "Unsure how to handle multiple bkbit categories:"
                        f"{json.dumps(record, indent=2)}"
                    )
                category = category[0]
                if category == 'bican:GenomeAnnotation':
                    candidate_name = record['id'].split(':')[-1]
                    candidate_name = candidate_name.replace(
                        'annotation-', ''
                    )
                    if citation_name is None:
                        citation_name = candidate_name
                    else:
                        if citation_name != candidate_name:
                            raise RuntimeError(
                                "More than one GenomeAnnotation in "
                                f"{bkbit_path}: ({citation_name} and "
                                f"{candidate_name}). Unsure how to proceed."
                            )
                if category in metadata_dict:
                    raise RuntimeError(
                        f"multiple entries for category {category}"
                    )
                metadata_dict[category] = record
                if 'authority' in record:
                    idx = metadata_utils._get_authority(
                        cursor=cursor,
                        name=record['authority'],
                        strict=True
                    )['idx']
                    if authority_idx is None:
                        authority_idx = idx
                    else:
                        if authority_idx != idx:
                            raise RuntimeError(
                                "Conflicting entries for authority"
                            )

    return (
        metadata_dict,
        raw_gene_annotation_values,
        citation_name
    )
