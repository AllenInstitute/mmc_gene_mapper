import pytest

from mmc_gene_mapper.utils.file_utils import (
    clean_up)


@pytest.fixture(scope='session')
def tmp_dir_fixture(
        tmp_path_factory):
    result = tmp_path_factory.mktemp('cell_type_mapper_')
    yield result
    clean_up(result)
