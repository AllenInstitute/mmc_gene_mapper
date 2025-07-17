import re


def int_from_identifier(identifier):
    """
    Take a gene identifier (a string) and return
    the integer part from the end as an integer
    """
    pattern = re.compile('[0-9]+$')
    result = pattern.findall(identifier)
    if len(result) != 1:
        raise MalformedGeneIdentifierError(
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
        assn = None
        if len(ival) > 0:
            ival = ival[0]
            ens_prefix = ens_prefix_pattern.match(gene_id)
            ncbi_prefix = ncbi_prefix_pattern.match(gene_id)
            if ens_prefix is not None and len(ival) > 0:
                if gene_id == f'{ens_prefix.group(0)}{ival}':
                    assn = 'ENSEMBL'
            elif ncbi_prefix is not None and len(ival) > 0:
                if gene_id == f'{ncbi_prefix.group(0)}{ival}':
                    assn = 'NCBI'
        if assn is None:
            assn = 'symbol'
        result.append(assn)
    return result


class MalformedGeneIdentifierError(Exception):
    pass
