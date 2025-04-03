"""
Define the class to actually do gene mapping
"""
import pathlib
import tempfile

import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.download.download_manager as download_manager
import mmc_gene_mapper.mapper.mapper_utils as mapper_utils


class MMCGeneMapper(object):

    def __init__(
           self,
           db_path,
           local_dir,
           data_file_spec=None,
           clobber=False,
           force_download=False):
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
        """

        dst_dir = pathlib.Path(
            tempfile.mkdtemp(
                dir=local_dir,
                prefix=(
                    f'mmc_gene_mapper_downloads_{timestamp.get_timestamp()}_'
                )
            )
        )

        self.download_mgr = download_manager.DownloadManager(
            dst_dir=dst_dir
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

            tmp_dir = tempfile.mkdtemp(
                dir=dst_dir,
                prefix='scratch_'
            )

            try:
                mapper_utils.create_mapper_database(
                    db_path=self.db_path,
                    download_manager=self.download_mgr,
                    data_file_spec=data_file_spec,
                    tmp_dir=tmp_dir,
                    force_download=force_download
                )
            finally:
                file_utils.clean_up(tmp_dir)

        contents = [n for n in dst_dir.iterdir()]
        if len(contents) == 0:
            file_utils.clean_up(dst_dir)
