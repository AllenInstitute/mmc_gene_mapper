"""
tests/with_data contains tests that run off of a cartoon data model

The structure of this cartoon data may be a little too complicated.
I am leaving it here because, but the time I realized this, there
wasn't a lot of time to refactor it and fix all the test cases
built on top of it. I will create a more systematic set of test
data in the arbitrary_conversion/ test module. It will be
tech debt to migrate all the tests in with_data/ over to using
that set of test data.

Apologies.
"""
import pytest

import hashlib
import itertools
import json
import pandas as pd
import pathlib
import tarfile
import tempfile
import unittest.mock


import mmc_gene_mapper.utils.file_utils as file_utils
import mmc_gene_mapper.mapper.mapper as mapper


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
def ncbi_ortholog_group_fixture():
    """
    List of lists of ints.
    Each list is a group of orthologous genes
    """
    return [
        [0, 14, 21],
        [1, 13, 27],
        [4, 12, 23],
        [7, 15],
        [6, 24]
    ]


@pytest.fixture(scope='session')
def alt_ortholog_group_fixture():
    """
    Alternative grouping of orthologous genes
    so that we can test querying orthologs from different
    citations.
    """
    return [
        [0, 11, 25],
        [2, 26, 13],
        [7, 12],
    ]


@pytest.fixture(scope='session')
def alt_ortholog_file_fixture(
        alt_ortholog_group_fixture,
        tmp_dir_fixture):
    """
    Write alternative ortholog groupings to
    a CSV file. Return the path to that file.
    """

    dst_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix="alternative_orthologs_",
        suffix=".csv"
    )

    data = []
    ct = 0
    for row in alt_ortholog_group_fixture:
        for g1 in row:
            data.append(
                {'ncbi_id': g1,
                 'garbage': ct,
                 'ortholog_id': row[-1],
                 }
            )
            ct += 1

    pd.DataFrame(data).to_csv(dst_path, index=False)
    return dst_path


@pytest.fixture(scope='session')
def ncbi_data_package_fixture(
        species_id_fixture,
        ncbi_ortholog_group_fixture):
    """
    Simulate data from which to generate a simulated NCBI-like test data
    package. Return a list of gene dicts each like
    {
     'species': int,
     'gene_id': int,
     'symbol': str,
     'ensembl_identifier': str
     'orthologs': list_of_ints
    }
    """
    species_list = [
        species_id_fixture['human'],
        species_id_fixture['jabberwock'],
        species_id_fixture['mouse']
    ]

    ortholog_lookup = dict()
    for row in ncbi_ortholog_group_fixture:
        for pair in itertools.combinations(row, 2):
            p0 = pair[0]
            p1 = pair[1]
            if p0 not in ortholog_lookup:
                ortholog_lookup[p0] = []
            ortholog_lookup[p0].append(p1)
            if p1 not in ortholog_lookup:
                ortholog_lookup[p1] = []
            ortholog_lookup[p1].append(p0)

    result = []
    for gene_id in range(38):
        # there will be mouse genes that are not listed in the bkbit data
        species_id = species_list[min(gene_id//10, 2)]

        symbol = f'symbol:{gene_id}'

        if gene_id == 5:
            symbol = "symbol:7"

        if gene_id % 2 == 1:
            ensembl = f'ENSX{2*gene_id}'
            ensembl_symbol = f'symbol:{2*(gene_id-1)}'
        else:
            ensembl = None
            ensembl_symbol = None

        if gene_id in ortholog_lookup:
            orthologs = ortholog_lookup[gene_id]
        else:
            orthologs = []
        result.append(
            {'species': species_id,
             'gene_id': gene_id,
             'symbol': symbol,
             'ensembl_identifier': ensembl,
             'ensembl_symbol': ensembl_symbol,
             'orthologs': orthologs}
        )

    return result


@pytest.fixture(scope='session')
def ncbi_file_package_fixture(
        ncbi_data_package_fixture,
        ncbi_ortholog_group_fixture,
        species_id_fixture,
        tmp_dir_fixture):
    """
    Return paths to files meant to look like an NCBI file package
    """
    tmp_dir = pathlib.Path(
        tempfile.mkdtemp(
            dir=tmp_dir_fixture,
            prefix='ncbi_package_'
        )
    )
    gene_info_path = tmp_dir / 'gene_info.gz'
    gene_2_ensembl_path = tmp_dir / 'gene2ensembl.gz'
    gene_ortholog_path = tmp_dir / 'gene_orthologs.gz'

    gene_info_data = [
        {'#tax_id': g['species'],
         'GeneID': g['gene_id'],
         'Symbol': g['symbol']}
        for g in ncbi_data_package_fixture
    ]
    pd.DataFrame(gene_info_data).to_csv(
        gene_info_path,
        sep='\t',
        compression='gzip',
        index=False
    )

    gene_2_ensembl_data = [
        {'#tax_id': g['species'],
         'GeneID': g['gene_id'],
         'Ensembl_gene_identifier': g['ensembl_identifier']}
        for g in ncbi_data_package_fixture
        if g['ensembl_identifier'] is not None
    ]

    # Add a gene whose ENSEMBL ID breaks the pattern
    # this should be gracefully skipped
    gene_2_ensembl_data.append(
        {'#tax_id': 10090,
         'GeneID': '1',
         'Ensembl_gene_identifier': 'ENS6666abcde'}
    )

    # simulate one ENSEMBL gene being equivalent
    # to more than one an NCBI gene
    gene_2_ensembl_data.append(
        {'#tax_id': species_id_fixture['human'],
         'GeneID': 5,
         'Ensembl_gene_identifier': 'ENSX14'}
    )

    pd.DataFrame(gene_2_ensembl_data).to_csv(
        gene_2_ensembl_path,
        sep='\t',
        compression='gzip',
        index=False
    )

    ortholog_data = []
    for row in ncbi_ortholog_group_fixture:
        for ii in range(1, len(row), 1):
            ortholog_data.append(
                {'GeneID': row[ii],
                 'Other_GeneID': row[ii-1],
                 'relationship': 'Ortholog'}
            )
    pd.DataFrame(ortholog_data).to_csv(
        gene_ortholog_path,
        sep='\t',
        compression='gzip',
        index=False
    )

    return {
        'gene_info': str(gene_info_path),
        'gene_2_ensembl': str(gene_2_ensembl_path),
        'gene_orthologs': str(gene_ortholog_path)
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


@pytest.fixture(scope='session')
def bkbit_data_fixture0(
        ncbi_data_package_fixture,
        tmp_dir_fixture):
    """
    Write out a simulated bkbit file. Return the path to
    that file.
    """
    json_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='test_bkbit_',
        suffix='.jsonld'
    )

    graph = []

    taxon = {
        "id": "NCBITaxon:999",
        "category": ["biolink:OrganismTaxon"],
        "name": "jabberwock",
        "full_name": "jabberwock"
    }

    assembly = {
        "id": "jabberwockAssembly",
        "category": ["bican:GenomeAssembly"],
        "name": "J001"
    }

    genome_annotation = {
        "id": "J001-2025",
        "category": ["bican:GenomeAnnotation"],
        "version": "0",
        "authority": "ENSEMBL"
    }

    ct = 0
    graph = [taxon, assembly, genome_annotation]
    for gene in ncbi_data_package_fixture:
        if gene['species'] != 999:
            continue
        if gene['ensembl_identifier'] is None:
            continue

        # give genes different symbols in ENSEMBL
        # than NCBI
        symbol = gene['ensembl_symbol']
        if ct % 2 == 1:
            name = f"name:{symbol.split(':')[-1]}"
        else:
            name = symbol
        ct += 1

        print(gene['ensembl_identifier'], symbol)

        gene = {
            "category": ["bican:GeneAnnotation"],
            "source_id": gene['ensembl_identifier'],
            "symbol": symbol,
            "name": name,
            "in_taxon_label": "jabberwock"
        }
        graph.append(gene)

    data = {'@graph': graph}
    with open(json_path, 'w') as dst:
        dst.write(json.dumps(data))
    return json_path


@pytest.fixture(scope='session')
def bkbit_data_fixture1(
        ncbi_data_package_fixture,
        tmp_dir_fixture):
    """
    Write out a simulated bkbit file for mouse. Return the path to
    that file.
    """
    json_path = file_utils.mkstemp_clean(
        dir=tmp_dir_fixture,
        prefix='test_bkbit_',
        suffix='.jsonld'
    )

    graph = []

    taxon = {
        "id": "NCBITaxon:10090",
        "category": ["biolink:OrganismTaxon"],
        "name": "mouse",
        "full_name": "Mus musculus"
    }

    assembly = {
        "id": "mouseAssembly",
        "category": ["bican:GenomeAssembly"],
        "name": "M001"
    }

    genome_annotation = {
        "id": "M001-2025",
        "category": ["bican:GenomeAnnotation"],
        "version": "0",
        "authority": "ENSEMBL"
    }

    ct = 0
    graph = [taxon, assembly, genome_annotation]
    for gene in ncbi_data_package_fixture:
        if gene['species'] != 10090:
            continue
        if gene['ensembl_identifier'] is None:
            continue
        idx = gene['gene_id']
        if idx >= 30:
            continue

        # give genes different symbols in ENSEMBL
        # than NCBI
        symbol = gene['ensembl_symbol']
        if ct % 2 == 1:
            name = f"name:{symbol.split(':')[-1]}"
        else:
            name = symbol
        ct += 1

        print(gene['ensembl_identifier'], symbol)

        gene = {
            "category": ["bican:GeneAnnotation"],
            "source_id": gene['ensembl_identifier'],
            "symbol": symbol,
            "name": name,
            "in_taxon_label": "Mus musculus"
        }
        graph.append(gene)

    data = {'@graph': graph}
    with open(json_path, 'w') as dst:
        dst.write(json.dumps(data))
    return json_path


@pytest.fixture(scope='session')
def dummy_download_mgr_fixture(
        species_file_fixture,
        ncbi_file_package_fixture):
    """
    define dummy downloan manager class
    """
    class DummyDownloadManager(object):
        def __init__(self, *args, **kwargs):
            pass

        def get_file(
                self,
                host,
                src_path,
                force_download):
            if src_path.endswith('new_taxdump.tar.gz'):
                return {
                    'local_path': species_file_fixture['tar']
                }
            elif src_path.endswith('new_taxdump.tar.gz.md5'):
                return {
                    'local_path': species_file_fixture['hash']
                }
            elif src_path.endswith('gene_info.gz'):
                return {
                    'local_path': ncbi_file_package_fixture['gene_info']
                }
            elif src_path.endswith('gene2ensembl.gz'):
                return {
                    'local_path': ncbi_file_package_fixture['gene_2_ensembl']
                }
            elif src_path.endswith('gene_orthologs.gz'):
                return {
                    'local_path': ncbi_file_package_fixture['gene_orthologs']
                }
            else:
                raise RuntimeError(
                    "MockDownloadManager cannot handle src_path "
                    f"{src_path}"
                )

    return DummyDownloadManager


@pytest.fixture(scope='session')
def mapper_fixture(
        ncbi_file_package_fixture,
        dummy_download_mgr_fixture,
        bkbit_data_fixture0,
        bkbit_data_fixture1,
        alt_ortholog_file_fixture,
        tmp_dir_fixture):
    """
    Return an instantiation of the MMCGeneMapper class
    based on our simulated NCBI file package
    """
    tmp_dir = tempfile.mkdtemp(dir=tmp_dir_fixture)
    db_path = file_utils.mkstemp_clean(
        dir=tmp_dir,
        prefix='ncbi_gene_mapper_',
        suffix='.db',
        delete=True
    )

    data_file_spec = [
        {"type": "bkbit",
         "path": bkbit_data_fixture0
         },
        {"type": "bkbit",
         "path": bkbit_data_fixture1},
        {"type": "hmba_orthologs",
         "path": alt_ortholog_file_fixture,
         "name": "alternative_orthologs"
         }
    ]

    to_replace = "mmc_gene_mapper.download.download_manager.DownloadManager"
    with unittest.mock.patch(to_replace, dummy_download_mgr_fixture):
        gene_mapper = mapper.MMCGeneMapper.create_mapper(
            db_path=db_path,
            local_dir=tmp_dir_fixture,
            data_file_spec=data_file_spec,
            clobber=False,
            force_download=False,
            suppress_download_stdout=True
        )

    return gene_mapper
