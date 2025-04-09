import pytest
import tempfile

import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.create_db.utils as db_utils


def test_check_existence(tmp_dir_fixture):

    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_fixture)
    tmp_file = file_utils.mkstemp_clean(dir=tmp_dir_fixture)
    assert db_utils.check_existence(tmp_file)
    assert not db_utils.check_existence('absolute_garbage.txt')
    with pytest.raises(file_utils.NotAFileError, match='is not a file'):
        db_utils.check_existence(tmp_dir)
