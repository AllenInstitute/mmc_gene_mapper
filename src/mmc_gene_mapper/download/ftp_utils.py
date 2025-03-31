"""
Define generic FTP download utils
"""
import ftplib
import json
import pathlib

import mmc_gene_mapper.utils.timestamp as timestamp



def download_files_from_ftp(
        ftp_host,
        file_dst_mapping,
        metadata_dst):
    """
    Parameters
    ----------
    ftp_host:
        string. The name of the FTP host from which to download
        the files
    file_dst_mapping:
        a dict mapping the names of the files on the FTP host to
        their destinations on the local drive
    metadata_dst:
        path to JSON file where metadata characterizing this
        download will be written
    """
    for dst in file_dst_mapping.values():
        dst = pathlib.Path(dst)
        if dst.exists():
            raise RuntimeError(
                f"{dst} already exists"
            )

    metadata = {
        'host': ftp_host,
        'files': list(file_dst_mapping.keys()),
        'downloaded_on': timestamp.get_timestamp()
    }
    metadata_str = json.dumps(metadata, indent=2)

    host = ftplib.FTP(ftp_host)
    host.login()
    for src in file_dst_mapping:
        print(f'=======DOWNLOADING {src}=======')
        with open(file_dst_mapping[src], 'wb') as dst:
            host.retrbinary(f'RETR {src}', dst.write)

    with open(metadata_dst, 'w') as dst:
        dst.write(metadata_str)
