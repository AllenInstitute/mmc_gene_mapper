"""
Utilities to help with managing temporary files
"""

import os
import pathlib
import tempfile


def assert_is_file(file_path):
    file_path = pathlib.Path(file_path)
    if not file_path.is_file():
        raise RuntimeError(
            f"{file_path} is not a file"
        )

def clean_up(target_path):
    """
    Recursively clean up and remove the directory at
    target_path
    """
    if target_path is None:
        return
    target_path = pathlib.Path(target_path)
    if target_path.is_file():
        target_path.unlink()
    elif target_path.is_dir():
        for sub_path in target_path.iterdir():
            clean_up(sub_path)
        target_path.rmdir()


def mkstemp_clean(
        dir=None,
        prefix=None,
        suffix=None,
        delete=False) -> str:
    """
    A thin wrapper around tempfile mkstemp that automatically
    closes the file descripter returned by mkstemp.

    Parameters
    ----------
    dir: Optional[Union[pathlib.Path, str]]
        The directory where the tempfile is created

    prefix: Optional[str]
        The prefix of the tempfile's name

    suffix: Optional[str]
        The suffix of the tempfile's name

    delete:
        if True, delete the file
        (dangerous, as could interfere with tempfile.mkstemp's
        ability to create unique file names; should only be used
        during testing where it is important that the tempfile
        doesn't actually exist)

    Returns
    -------
    file_path: str
        Path to a valid temporary file

    Notes
    -----
    Because this calls tempfile mkstemp, the file will be created,
    though it will be empty. This wrapper is needed because
    mkstemp automatically returns an open file descriptor, which was
    been causing some of our unit tests to overwhelm the OS's limit
    on the number of open files.
    """
    (descriptor,
     file_path) = tempfile.mkstemp(
                     dir=dir,
                     prefix=prefix,
                     suffix=suffix)

    os.close(descriptor)
    if delete:
        os.unlink(file_path)
    return file_path
