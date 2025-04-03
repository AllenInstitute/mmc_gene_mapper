import re


def int_from_identifier(identifier):
    pattern = re.compile('[0-9]+')
    result = pattern.findall(identifier)
    if len(result) != 1:
        raise ValueError(
            f'could not get one integer from {identifier}'
        )
    return int(result[0])


def detect_gene_identifiers(
        gene_id_list):
    """
    Take a list of gene identifiers.

    If we think they are identifiers return the
    authority we think they belong to.

    If we think they are symbols, return None

    Must make the same choice for all identifiers
    """
    try:
        [int_from_identifier(ii) for ii in gene_id_list]
    except ValueError:
        return None

    ncbi_start = [
        ii.startswith('NCBI') for ii in gene_id_list
    ]
    if all(ncbi_start):
        return 'NCBI'
    ens_start = [
        ii.startswith('ENS') for ii in gene_id_list
    ]
    if all(ens_start):
        return 'ENSEMBL'
    return None
