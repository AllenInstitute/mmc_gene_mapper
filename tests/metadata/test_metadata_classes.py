import pytest

import mmc_gene_mapper.metadata.classes as metadata_classes


def test_authority_errors():
    """
    Make sure invalid authority raises an error
    """
    msg = "authority = 'gla' is not a valid"
    with pytest.raises(ValueError, match=msg):
        metadata_classes.Authority('gla')
