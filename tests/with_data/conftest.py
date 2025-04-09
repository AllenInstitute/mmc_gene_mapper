"""
tests/with_data contains tests that run off of a cartoon data model
"""
import pytest

import hashlib
import pathlib
import tarfile
import tempfile


@pytest.fixture(scope='session')
def species_id_fixture():
    """
    Return a dict mapping species names to taxon_ids
    """
    return {
        'mouse': 10090,
        'Mus musculus': 10090,
        'human': 9606,
        'Homo sapiens': 9606,
        'jabberwock': 999,
        'elf': 998,
        'fey': 998,
        'goblin': 997,
        'orc': 997
    }


@pytest.fixture(scope='session')
def species_file_fixture(
        species_id_fixture,
        tmp_dir_fixture):
    """
    Create a txv simulating names.dmp. Put it into a tarfile.
    Create the text file with that tarfile's md5 hash.
    Return the path to the tarfile and the text file.
    This is meant to simulate downloading the taxons from NCBI.
    """
    tmp_dir = pathlib.Path(
        tempfile.mkdtemp(
            dir=tmp_dir_fixture,
            prefix='species_simulation_'
        )
    )

    name_path = tmp_dir / 'names.dmp'
    with open(name_path, 'w') as dst:
        for species_name in species_id_fixture:
            species_id = species_id_fixture[species_name]
            dst.write(f'{species_id}\t|\t{species_name}\t|\t\n')

    tar_path = tmp_dir / 'new_taxdump.tar.gz'
    with tarfile.open(tar_path, 'w:gz') as dst:
        with open(name_path, 'rb') as src:
            t_info = dst.gettarinfo(fileobj=src)
            t_info.name = name_path.name
            dst.addfile(
                tarinfo=t_info,
                fileobj=src
            )

    hasher = hashlib.md5()
    with open(tar_path, 'rb') as src:
        hasher.update(src.read())
    hash_path = tmp_dir / 'new_taxdump.tar.gz.md5'
    with open(hash_path, 'w') as dst:
        dst.write(hasher.hexdigest())

    return {"tar": tar_path, "hash": hash_path}
