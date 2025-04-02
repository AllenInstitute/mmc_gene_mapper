import re


def int_from_identifier(identifier):
    pattern = re.compile('[0-9]+')
    result = pattern.findall(identifier)
    if len(result) != 1:
        raise ValueError(
            f'could not get one integer from {identifier}'
        )
    return int(result[0])
