"""
Define generic FTP download utils
"""
import json
import pathlib
import subprocess
import time

import mmc_gene_mapper.utils.timestamp as timestamp
import mmc_gene_mapper.utils.file_utils as file_utils


def download_file(
        host,
        src_path,
        dst_path,
        clobber=False,
        suppress_stdout=False):
    """
    Parameters
    ----------
    host:
        string. The name of the FTP host from which to download
        the files
    src_path:
        path within the host of the file being downloaded
    dst_path:
        path on local machine where we will be saving the file
    clobber:
        whether or not to overwrite existing file
    suppress_stdout:
        if True, pipe stdout and stderr to /dev/NULL

    Returns
    -------
    metadata describing downloaded file
    """

    dst_path = pathlib.Path(dst_path)
    if dst_path.exists():
        if clobber:
            file_utils.clean_up(dst_path)
        else:
            raise RuntimeError(
                f"{dst_path} already exists"
            )

    metadata = {
        'host': host,
        'src_path': src_path,
        'downloaded_on': timestamp.get_timestamp(),
        'local_path': str(dst_path.resolve().absolute())
    }
    json.dumps(metadata, indent=2)

    print(f'=======DOWNLOADING {src_path}=======')
    src_url = (
        f'https://{host}/{src_path}'
    )
    dst_url = (
        f'{str(dst_path.resolve().absolute())}'
    )
    process_args = [
        "wget",
        src_url,
        "-O",
        dst_url
    ]

    if suppress_stdout:
        process_args.append("-q")

    t0 = time.time()
    process = subprocess.Popen(
        args=process_args)
    return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(
            f"subprocess {process_args} returned code {return_code}"
        )
    dur = time.time()-t0
    print(f'=======DOWNLOADED {src_path} in {dur:.2e} seconds=======')
    return metadata
