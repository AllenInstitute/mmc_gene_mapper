"""
Define generic FTP download utils
"""
import json
import pathlib
import subprocess

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

    for src in file_dst_mapping:
        print(f'=======DOWNLOADING {src}=======')
        src_url = (
            f'https://{ftp_host}/{src}'
        )
        dst = pathlib.Path(file_dst_mapping[src])
        dst_url = (
            f'{str(dst.resolve().absolute())}'
        )
        process_args=[
            "wget",
            src_url,
            "-O",
            dst_url
        ]
        process = subprocess.Popen(args=process_args)
        return_code = process.wait()
        if return_code != 0:
            raise RuntimeError(
                f"subprocess {args} returned code {return_code}"
            )

    with open(metadata_dst, 'w') as dst:
        dst.write(metadata_str)
