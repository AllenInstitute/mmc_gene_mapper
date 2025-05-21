"""
Define functions and classes to manage downloading of
files
"""

import pathlib

import mmc_gene_mapper
import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.download.download_utils as download_utils
import mmc_gene_mapper.download.download_manager_utils as mgr_utils


class DownloadManager(object):

    def __init__(
            self,
            dst_dir,
            suppress_stdout=False):

        self.suppress_stdout = suppress_stdout

        dst_dir = pathlib.Path(dst_dir)
        if not dst_dir.is_dir():
            try:
                dst_dir.mkdir(parents=True)
            except Exception:
                raise RuntimeError(
                    f"{dst_dir} is not a dir; "
                    "cannot download data there"
                )
        self.dst_dir = dst_dir

        data_dir = pathlib.Path(
            mmc_gene_mapper.__file__
        ).parent / "download/db_files"
        if not data_dir.is_dir():
            raise RuntimeError(
                f"{data_dir} is not a dir"
            )
        self.db_path = data_dir / 'download_manager.db'
        if self.db_path.exists():
            if not self.db_path.is_file():
                raise RuntimeError(
                    f"{self.db_path} is not a file"
                )
        else:
            mgr_utils.create_download_db(self.db_path)

    def remove_file(
            self,
            host,
            src_path):

        mgr_utils.remove_record(
            db_path=self.db_path,
            host=host,
            src_path=src_path
        )

    def get_file(
            self,
            host,
            src_path,
            force_download=True):

        pre_existing = mgr_utils.get_record(
            db_path=self.db_path,
            host=host,
            src_path=src_path
        )
        if len(pre_existing) == 0:
            do_download = True
        elif len(pre_existing) == 1 and not force_download:
            record = pre_existing[0]
            local_path = pathlib.Path(
                record['local_path']
            )
            if not local_path.is_file():
                do_download = True
            else:
                actual_hash = file_utils.hash_from_path(
                    local_path
                )
                if actual_hash != record['hash']:
                    do_download = True
                else:
                    return record
        else:
            do_download = True

        if do_download:
            if len(pre_existing) > 0:
                mgr_utils.remove_record(
                    db_path=self.db_path,
                    host=host,
                    src_path=src_path
                )

            fname = pathlib.Path(src_path).name
            suffix = pathlib.Path(src_path).suffix
            salt = 0
            dst_path = None
            while True:
                if suffix == '':
                    dst_path = self.dst_dir / f'{fname}.{salt}'
                else:
                    dst_path = self.dst_dir / fname.replace(
                        suffix, f'.{salt}{suffix}'
                    )
                if not dst_path.exists():
                    break
                salt += 1
            download_metadata = download_utils.download_file(
                host=host,
                src_path=src_path,
                dst_path=dst_path,
                clobber=False,
                suppress_stdout=self.suppress_stdout
            )
            mgr_utils.insert_record(
                db_path=self.db_path,
                host=download_metadata['host'],
                src_path=download_metadata['src_path'],
                local_path=download_metadata['local_path']
            )

        record = mgr_utils.get_record(
            db_path=self.db_path,
            host=host,
            src_path=src_path
        )
        if len(record) != 1:
            raise RuntimeError(
                "Even after forced download, more than one record for "
                f"host: {host}, src_path: {src_path}"
            )
        return record[0]
