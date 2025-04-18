"""
Define the class to actually do gene mapping
"""
import json
import pathlib
import shutil
import sqlite3
import tempfile
import time

import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.download.download_manager as download_manager
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils
import mmc_gene_mapper.query_db.query as query_utils


class MMCGeneMapper(object):

    def __init__(
           self,
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
            self._initialize(
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

    def identifiers_from_symbols(
            self,
            gene_symbol_list,
            species_name,
            authority_name):
        """
        Find the mapping that converts gene symbols into
        gene identifiers

        Parameters
        ----------
        gene_symbol_list:
            list of gene symbols
        species_name:
            name of the species we are working with
        authority_name:
            name of the authority in whose identifiers
            we want the genes listed

        Returns
        -------
        A dict
            {
              "metadata": {
                  a dict describing the citation according
                  to which these symbols map to these
                  identifiers
              },
              "mapping": {
                  a dict keyed on input symbols. Each symbol
                  maps to the list of all gene identifiers
                  that are associated with that symbol according
                  to the source described in "metadata"
              }
            }
        """
        species_taxon = query_utils.get_species_taxon(
            db_path=self.db_path,
            species_name=species_name,
            strict=True)

        result = query_utils.translate_gene_identifiers(
            db_path=self.db_path,
            src_column="symbol",
            dst_column="identifier",
            src_list=gene_symbol_list,
            authority_name=authority_name,
            species_taxon=species_taxon,
            chunk_size=500
        )
        return result

    def equivalent_genes(
            self,
            input_authority,
            output_authority,
            gene_list,
            species_name,
            citation_name):
        """
        Return a mapping between gene identifiers from
        different authorities (NCBI vs ENSEMBL)

        Parameters
        ----------
        input_authority:
            a str; the name of the authority in which the input
            genes are identified
        output_authority:
            a str; the name of the authority you want to convert
            the identifiers to
        gene_list:
            list of gene identifiers (in input_authority) to be
            mapped
        species_name:
            name of species we are working with
        citation_name:
            name of citation to use to assess gene eqivalence

        Returns
        -------
        A dict
            {
              "metadata": {
                  a dict describing the citation according
                  to which these symbols map to these
                  identifiers
              },
              "mapping": {
                  a dict keyed on input genes. Each gene
                  maps to the list of all genes
                  that are considered equivalent to that
                  gene according to the source described
                  in "metadata"
              }
            }
        """
        return query_utils.get_equivalent_genes_from_identifiers(
            db_path=self.db_path,
            input_authority_name=input_authority,
            output_authority_name=output_authority,
            input_gene_list=gene_list,
            species_name=species_name,
            citation_name=citation_name,
            chunk_size=500
        )

    def ortholog_genes(
            self,
            authority,
            src_species_name,
            dst_species_name,
            gene_list,
            citation_name):
        """
        Return a mapping between gene identifiers from
        different species

        Parameters
        ----------
        authority:
            a str; the name of the authority (ENSEMBL, NCBI etc.)
            we are working in
        src_species_name:
            a str; the name of the species we are starting from
        dst_species_name:
            as str; the name of the species we are mapping to
        gene_list:
            list of gene identifiers (in src_species) to be
            mapped
        citation_name:
            name of citation to use to assess gene eqivalence

        Returns
        -------
        A dict
            {
              "metadata": {
                  a dict describing the citation according
                  to which these symbols map to these
                  identifiers
              },
              "mapping": {
                  a dict keyed on input genes. Each gene
                  maps to the list of all genes
                  that are considered orthologs to that
                  gene according to the source described
                  in "metadata"
              }
            }
        """
        return query_utils.get_ortholog_genes_from_identifiers(
            db_path=self.db_path,
            authority_name=authority,
            src_species_name=src_species_name,
            dst_species_name=dst_species_name,
            src_gene_list=gene_list,
            citation_name=citation_name,
            chunk_size=500
        )

    def _initialize(
            self,
            dst_dir,
            tmp_dir,
            db_path,
            data_file_spec,
            clobber,
            force_download,
            suppress_download_stdout):

        self.download_mgr = download_manager.DownloadManager(
            dst_dir=dst_dir,
            suppress_stdout=suppress_download_stdout
        )

        self.db_path = pathlib.Path(db_path)
        if self.db_path.exists():
            if not self.db_path.is_file():
                raise RuntimeError(
                    f"db_path {self.db_path} is not a file"
                )
            if clobber:
                self.db_path.unlink()

        if not self.db_path.is_file():

            t0 = time.time()
            print('=======CREATING DB FILE=======')
            tmp_db_path = file_utils.mkstemp_clean(
                dir=tmp_dir,
                suffix='.db',
                delete=True
            )

            mapper_utils.create_mapper_database(
                db_path=tmp_db_path,
                download_manager=self.download_mgr,
                data_file_spec=data_file_spec,
                tmp_dir=tmp_dir,
                force_download=force_download
            )

            mapper_utils.create_bibliography_table(
                tmp_db_path
            )

            print(f'=======COPYING TMP FILE TO {self.db_path}=======')
            shutil.move(
                src=tmp_db_path,
                dst=self.db_path
            )
            dur = (time.time()-t0)/60.0
            print(f'=======DB CREATION TOOK {dur:.2e} MINUTES=======')
