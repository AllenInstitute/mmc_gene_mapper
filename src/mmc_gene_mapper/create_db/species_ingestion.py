"""
Functions for populating the species table
"""
import hashlib
import pathlib
import sqlite3
import tarfile
import tempfile

import mmc_gene_mapper.utils.file_utils as file_utils


def ingest_species_data(
        db_path,
        download_manager,
        force_download=False,
        tmp_dir=None):

    tmp_dir = pathlib.Path(
        tempfile.mkdtemp(dir=tmp_dir)
    )

    if not tmp_dir.is_dir():
        raise RuntimeError(
            f"{tmp_dir} is not a dir"
        )

    try:
        _ingest_species_data(
            db_path=db_path,
            download_manager=download_manager,
            force_download=force_download,
            tmp_dir=tmp_dir
        )
    finally:
        file_utils.clean_up(tmp_dir)


def _ingest_species_data(
        db_path,
        download_manager,
        force_download,
        tmp_dir):

    db_path = pathlib.Path(db_path)
    if db_path.exists():
        file_utils.assert_is_file(db_path)

    host = 'ftp.ncbi.nih.gov'
    tar_record = download_manager.get_file(
        host=host,
        src_path='pub/taxonomy/new_taxdump/new_taxdump.tar.gz',
        force_download=force_download
    )

    md5_record = download_manager.get_file(
        host=host,
        src_path='pub/taxonomy/new_taxdump/new_taxdump.tar.gz.md5',
        force_download=force_download
    )

    tar_path = tar_record['local_path']
    md5_path = md5_record['local_path']

    md5 = hashlib.md5()
    with open(tar_path, 'rb') as src:
        while True:
            chunk = src.read(100000000)
            if len(chunk) == 0:
                break
            md5.update(chunk)

    with open(md5_path, 'r') as src:
        expected = src.readline().split()[0]

    if expected != md5.hexdigest():
        raise RuntimeError(
            "Hash for new_taxdump.tar.gz does not match"
        )

    name_path = tmp_dir / 'names.dmp'
    if name_path.exists():
        raise RuntimeError(
            f"{name_path} exists before untarring"
        )

    with tarfile.open(tar_path, mode='r') as src:
        src.extract('names.dmp', path=tmp_dir, filter='data')

    if not name_path.is_file():
        raise RuntimeError(
            f"{name_path} is not a file after untarring"
        )

    ingest_species_table(
        db_path=db_path,
        data_path=name_path)


def ingest_species_table(
        db_path,
        data_path):

    species = []
    with open(data_path, "r") as src:
        for line in src:
            params = [el.strip() for el in line.split('|')]
            species.append((int(params[0]), params[1]))

    table_name = "NCBI_species"
    index_name = "NCBI_species_idx"
    with sqlite3.connect(db_path) as conn:
        cursor = conn.cursor()
        cursor.execute(
            f"DROP INDEX IF EXISTS {index_name}"
        )
        cursor.execute(
            f"DROP TABLE IF EXISTS {table_name}"
        )

        cursor.execute(
            f"""
            CREATE TABLE {table_name} (
                id INTEGER,
                name STRING
            )
            """
        )

        cursor.executemany(
            f"""
            INSERT INTO {table_name} (
                id,
                name
            )
            VALUES (?, ?)
            """,
            species
        )

        cursor.execute(
            f"""
            CREATE INDEX {index_name} on {table_name} (name)
            """
        )
