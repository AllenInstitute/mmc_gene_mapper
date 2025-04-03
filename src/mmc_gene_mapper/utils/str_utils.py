import re


def int_from_identifier(identifier):
    pattern = re.compile('[0-9]+')
    result = pattern.findall(identifier)
    if len(result) != 1:
        raise ValueError(
            f'could not get one integer from {identifier}'
        )
    return int(result[0])


def characterize_gene_identifiers(
        gene_id_list):
    """
    Take a list of gene identifiers.

    Return a list indicating if genes are
        ['NCBI', 'ENSEMBL', or 'symbol']

    """
    int_pattern = re.compile('[0-9]+')
    ens_prefix_pattern = re.compile('ENS[A-Z]+')
    ncbi_prefix_pattern = re.compile('NCBI[A-Za-z]*[:]?')
    result = []
    for gene_id in gene_id_list:
        ival = int_pattern.findall(gene_id)
        if len(ival) > 0:
            ival = ival[0]
        else:
            ival = ''
        ens_prefix = ens_prefix_pattern.match(gene_id)
        ncbi_prefix = ncbi_prefix_pattern.match(gene_id)
        if ens_prefix is not None:
            if gene_id == f'{ens_prefix.group(0)}{ival}':
                result.append('ENSEMBL')
                continue
        elif ncbi_prefix is not None:
            if gene_id == f'{ncbi_prefix.group(0)}{ival}':
                result.append('NCBI')
                continue
        result.append('symbol')
    return result
