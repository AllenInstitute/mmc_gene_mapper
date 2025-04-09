import pytest

import hashlib

import mmc_gene_mapper.utils.file_utils as file_utils


def test_assert_is_file(tmp_dir_fixture):

    with pytest.raises(file_utils.NotAFileError):
        file_utils.assert_is_file("absolute_garbage.csv")

    actual_file = file_utils.mkstemp_clean(dir=tmp_dir_fixture)
    file_utils.assert_is_file(actual_file)


def test_hash_from_path(tmp_dir_fixture):

    dst_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        suffix='.txt'
    )
    with open(dst_path, 'w') as dst:
        dst.write('some text')

    expected = hashlib.md5()
    with open(dst_path, 'rb') as src:
        expected.update(src.read())

    actual = file_utils.hash_from_path(dst_path)
    assert actual == f'md5:{expected.hexdigest()}'
