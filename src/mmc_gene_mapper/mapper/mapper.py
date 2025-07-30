"""
Define the class to actually do gene mapping
"""
import json
import pathlib
import shutil
import sqlite3
import tempfile
import traceback
import time

import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.metadata.classes as metadata_classes
import mmc_gene_mapper.download.download_manager as download_manager
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils
import mmc_gene_mapper.query_db.query as query_utils
import mmc_gene_mapper.mapper.arbitrary_conversion as arbitrary_conversion


class MMCGeneMapper(object):

    def __init__(
            self,
            db_path):
        self.db_path = pathlib.Path(db_path)
        if not self.db_path.is_file():
            raise MalformedMapperDBError(
                f"db_path {self.db_path} is not a file"
            )
        try:
            with sqlite3.connect(db_path) as conn:
                cursor = conn.cursor()
                validity = cursor.execute(
                    "SELECT validity FROM mmc_gene_mapper_metadata"
                ).fetchall()
                if validity != [('TRUE',),]:
                    raise MalformedMapperDBError(
                        f"file at {db_path} "
                        "not marked as valid mmc_gene_mapper "
                        f"database; validity = {validity}"
                    )
        except Exception:
            err_msg = f"\n{traceback.format_exc()}\n"
            raise MalformedMapperDBError(
                f"An error occurred while validating file at {db_path}"
                f"{err_msg}"
            )

    @classmethod
    def create_mapper(
           cls,
           db_path,
           local_dir,
           data_file_spec=None,
           clobber=False,
           force_download=False,
           suppress_download_stdout=False):
        """
        Parameters
        ----------
        db_path:
            path to the database where gene mapping data will be stored
        local_dir:
            path to directory to which files can be downloaded when
            creating a database (if necessary)
        data_file_spec:
            list of dicts like
                {'type': ('bkbit', or 'hmba_ortholog')
                 'name': name_by_which_to_refer_to_citation,
                'path': path to file}
        clobber:
            if True, overwrite existing database
        force_download:
            if True, re-download data
        suppress_download_stdout:
            if True, suppress the stdout produced by calling
            wget with subprocess while downloading data
            (this is necessary when running in a notebook)
        """

        dst_dir = pathlib.Path(
            tempfile.mkdtemp(
                dir=local_dir,
                prefix=(
                    f'mmc_gene_mapper_downloads_{timestamp.get_timestamp()}_'
                )
            )
        )

        tmp_dir = pathlib.Path(
            tempfile.mkdtemp(
                dir=dst_dir,
                prefix='scratch_'
            )
        )
        try:
            _initialize_mapper(
                dst_dir=dst_dir,
                tmp_dir=tmp_dir,
                db_path=db_path,
                data_file_spec=data_file_spec,
                clobber=clobber,
                force_download=force_download,
                suppress_download_stdout=suppress_download_stdout
            )
        finally:

            file_utils.clean_up(tmp_dir)

            contents = [n for n in dst_dir.iterdir()]
            if len(contents) == 0:
                file_utils.clean_up(dst_dir)

        return cls(db_path=db_path)

    def get_all_species(self):
        """
        Return a list of all species names in the database.

        Note: these are not all the species for which there are
        genes, orthologs, etc. These are simply all of the species
        which the database can translate into a numerial value
        for cross referencing against the gene tables.
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            gene_names = cursor.execute(
                """
                SELECT
                    DISTINCT(name)
                FROM
                    NCBI_species
                """
            ).fetchall()
        result = sorted([str(row[0]) for row in gene_names])
        return result

    def get_all_citations(self):
        """
        Return a list of dicts representing all of the citations
        recorded in this database.

        Dicts will be of the form
            {"name": "NameOfCitation",
             "metadata": {dict representing metadata describing citation}
            }
        """

        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            raw = cursor.execute(
                """
                SELECT
                    name,
                    metadata
                FROM
                    citation
                """
            ).fetchall()

        return [
            {"name": row[0],
             "metadata": json.loads(row[1])}
            for row in raw
        ]

    def get_all_authorities(self):
        """
        Return a list of strings representing the names of all
        the authorities in this database
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            raw = cursor.execute(
                """
                SELECT
                    name
                FROM
                    authority
                """
            ).fetchall()

        return [
            row[0] for row in raw
        ]

    def map_genes(
            self,
            gene_list,
            dst_species,
            dst_authority,
            ortholog_citation='NCBI',
            log=None,
            invalid_mapping_prefix=None):
        """
        Perform an arbitrary mapping on a list of gene identifiers.

        Parameters
        ----------
        db_path:
            path to the database being queries
        gene_list:
            list of gene identifiers being mapped
        dst_species:
            name of species being mapped to
        dst_authority:
            name of authority being mapped to
        ortholog_citation:
            citation to use for ortholog mapping, if necessary
        log:
            a logger class that implements an info()
            function (probably the CommandLog from cell_type_mapper)
        invalid_mapping_prefix:
            an optional string. If not None, this will be prepended
            to the placeholder names of all unmappable genes
        """
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            dst_species = query_utils.get_species(
                cursor=cursor,
                species=dst_species
            )

        return arbitrary_conversion.arbitrary_mapping(
            db_path=self.db_path,
            gene_list=gene_list,
            dst_species=dst_species,
            dst_authority=metadata_classes.Authority(dst_authority),
            ortholog_citation=ortholog_citation,
            log=log,
            invalid_mapping_prefix=invalid_mapping_prefix
        )


def _initialize_mapper(
        dst_dir,
        tmp_dir,
        db_path,
        data_file_spec,
        clobber,
        force_download,
        suppress_download_stdout):

    download_mgr = download_manager.DownloadManager(
        dst_dir=dst_dir,
        suppress_stdout=suppress_download_stdout
    )

    db_path = pathlib.Path(db_path)
    if db_path.exists():
        if not db_path.is_file():
            raise RuntimeError(
                f"db_path {db_path} is not a file"
            )
        if clobber:
            db_path.unlink()
        else:
            raise RuntimeError(
                f"db_path {db_path} already exists; run with clobber=True "
                "to delete and re-create"
            )

    t0 = time.time()
    print('=======CREATING DB FILE=======')
    tmp_db_path = file_utils.mkstemp_clean(
        dir=tmp_dir,
        suffix='.db',
        delete=True
    )

    mapper_utils.create_mapper_database(
        db_path=tmp_db_path,
        download_manager=download_mgr,
        data_file_spec=data_file_spec,
        tmp_dir=tmp_dir,
        force_download=force_download
    )

    mapper_utils.create_bibliography_table(
       tmp_db_path
    )

    pre_metadata_hash = file_utils.hash_from_path(
        file_path=tmp_db_path
    )

    # mark this as a valid mmc_gene_mapper database
    with sqlite3.connect(tmp_db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            """
            CREATE TABLE mmc_gene_mapper_metadata (
                validity STR,
                timestamp STR,
                hash STR
            )
            """
        )
        cursor.execute(
            """
            INSERT INTO mmc_gene_mapper_metadata
                (validity,
                 timestamp,
                 hash)
             VALUES("TRUE", ?, ?)
            """,
            (timestamp.get_timestamp(), pre_metadata_hash)
        )

    print(f'=======COPYING TMP FILE TO {db_path}=======')
    shutil.move(
        src=tmp_db_path,
        dst=db_path
    )
    dur = (time.time()-t0)/60.0
    print(f'=======DB CREATION TOOK {dur:.2e} MINUTES=======')


class MalformedMapperDBError(Exception):
    pass
