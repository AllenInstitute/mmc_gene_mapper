import pytest

import pathlib
import tempfile

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.download.download_manager_utils as mgr_utils


def test_create_download_db(tmp_dir_fixture):
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_fixture)
    existing_file = file_utils.mkstemp_clean(
        dir=tmp_dir,
        prefix='existing_db_',
        suffix='.db',
        delete=False
    )

    with pytest.raises(ValueError, match="already exists"):
        mgr_utils.create_download_db(existing_file)

    valid_path = pathlib.Path(tmp_dir) / 'not_a_file.db'
    assert not valid_path.exists()
    mgr_utils.create_download_db(valid_path)
    assert valid_path.is_file()
