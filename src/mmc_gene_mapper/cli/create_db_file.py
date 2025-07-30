import argparse
import pathlib
import tempfile

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.ensembl_download.scraper as ensembl_scraper
import mmc_gene_mapper.mapper.mapper as mapper_module


def main():
    parser = argparse.ArgumentParser(
        "Create the sqlite db file used by mmc_gene_mapper to "
        "map genes across species and between authorities "
        "(i.e. 'NCBI' and 'ENSEMBL')"
    )
    parser.add_argument(
       "--db_path",
       type=str,
       default=None,
       help="Path to the file that will be written."
    )
    parser.add_argument(
        "--local_dir",
        type=str,
        default=None,
        help=(
            "Directory to which data files can be downloaded. "
            "Must be empty. Will be created if does not exist."
        )
    )
    parser.add_argument(
        "--ensembl_version",
        type=int,
        default=114,
        help="Version of ENSEMBL to scrape (default=114)."
    )
    parser.add_argument(
        "--suppress_stdout",
        default=False,
        action="store_true",
        help=(
            "Whether or not to suppress stdout messages "
            "produced by file download."
        )
    )
    parser.add_argument(
        "--clobber",
        default=False,
        action="store_true",
        help="Whether or not to overwrite existing file at db_path."
    )

    args = parser.parse_args()

    create_db_file(
        db_path=args.db_path,
        local_dir=args.local_dir,
        ensembl_version=args.ensembl_version,
        suppress_download_stdout=args.suppress_stdout,
        clobber=args.clobber
    )


def create_db_file(
        db_path,
        local_dir,
        ensembl_version,
        clobber,
        suppress_download_stdout):
    """
    Create sqlite db file for backing mmc_gene_mapper

    Parameters
    ----------
    db_path:
        path to db file that will be written
    local_dir:
        path to local directory to which files can be downloaded.
        Must be empty. Will be created if it does not exist.
    clobber:
        boolean indicating whether or not to overwrite file that
        already exist at db_path
    suppress_download_stdout:
        boolean. If True, suppress stdout produced by file downloads
    """

    local_dir = pathlib.Path(local_dir)
    if local_dir.exists():
        if not local_dir.is_dir():
            raise ValueError(
                f"local_dir {local_dir} is not a dir"
            )
        contents = [n for n in local_dir.iterdir()]
        if len(contents) != 0:
            raise ValueError(
                f"local_dir {local_dir} is not empty"
            )
    else:
        local_dir.mkdir(parents=True)

    db_path = pathlib.Path(db_path)
    if db_path.exists():
        if not db_path.is_file():
            raise ValueError(f"db_path {db_path} is not a file")
        else:
            if clobber:
                db_path.unlink()
            else:
                raise ValueError(
                    f"db_path {db_path} already exists; "
                    "run with --clobber to overwrite"
                )

    bican_files = [
        {'url': (
            'https://ftp.ensembl.org/pub/release-101/gff3/'
            'homo_sapiens/Homo_sapiens.GRCh38.101.gff3.gz'
         ),
         'assembly_id': 'GCF_000001405.39'},
        {'url': (
            'https://ftp.ensembl.org/pub/release-98/gff3/'
            'mus_musculus/Mus_musculus.GRCm38.98.gff3.gz'
         ),
         'assembly_id': 'GCF_000001635.26'}
    ]

    scratch_dir = pathlib.Path(
        tempfile.mkdtemp(dir=local_dir, prefix='scratch')
    )
    ensembl_download_dir = pathlib.Path(
        tempfile.mkdtemp(dir=scratch_dir, prefix='ensembl')
    )
    try:
        print("====SCRAPING ENSEMBL====")
        if ensemb_version > 0:
            ensembl_files_spec = ensembl_scraper.scrape_ensembl(
                default_ensembl_version=ensembl_version,
                dst_dir=ensembl_download_dir,
                failure_log_path=file_utils.mkstemp_clean(
                    dir=scratch_dir,
                    suffix='.txt'
                ),
                specific_files=bican_files,
                tmp_dir=scratch_dir,
                n_limit=None
            )
            print("====DONE SCRAPING ENSEMBL====")
        else:
            ensembl_files_spec = None

        mapper_module.MMCGeneMapper.create_mapper(
            db_path=db_path,
            local_dir=local_dir,
            data_file_spec=ensembl_files_spec,
            clobber=clobber,
            force_download=True,
            suppress_download_stdout=suppress_download_stdout
        )

        print(
            "SUCCESS; wrote {str(db_path)}"
        )

    finally:
        file_utils.clean_up(scratch_dir)


if __name__ == "__main__":
    main()
