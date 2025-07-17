"""
Define functions to scrape ENSEMBL ftp server for data
"""
import copy
import ftplib
import json
import pathlib
import tempfile
import time
import traceback
import warnings

import mmc_gene_mapper.utils.file_utils as file_utils

import bkbit.data_translators.genome_annotation_translator as gene_translator


def scrape_ensembl(
        default_ensembl_version,
        dst_dir,
        failure_log_path,
        specific_files=None,
        n_limit=None,
        tmp_dir=None,
        clobber=False):
    """
    Parameters
    ----------
    default_ensembl_version:
        an int. The default release version to download from ENSEMBL
    dst_dir:
        directory where jsonld files can be written
    failure_log_path:
        path to file where files that failed to process can be downloaded
    specific_files:
        optional list of specifications of specific versions of species
        genome annotations to download. Each entry is a dict like
        {'url': http://url/to/file.gff3.gz,
         'assembly_id': 'id fo genome assembly'}
    n_limit:
        optional it to limit the number of files actually parsed
        (for prototyping)
    tmp_dir:
        a directory where temporary files downloaded from ENSEMBL can be
        saved while they are needed.
    clobber:
        boolean. If true and dst_dir exists, delete it before proceeding

    Returns
    -------
    A list of dicts. The data_file_spec to be passed to MMCGeneMapper for the
    jsonld files produced

    Notes
    -----
    if dst_dir exists and is not empty, this function will fail unless
    you run with clobber=True, emptying and reconstituting dst_dir before doing
    any processing.
    """
    t0 = time.time()
    dst_dir = pathlib.Path(dst_dir)
    if dst_dir.exists():
        if not dst_dir.is_dir():
            raise ValueError(
                f"{dst_dir} is not a dir"
            )

        contents = [n for n in dst_dir.iterdir()]
        if len(contents) > 0:
            if clobber:
                file_utils.clean_up(dst_dir)
            else:
                raise ValueError(
                    f"dst_dir={dst_dir} exists and is not empty; "
                    "run with clobber=True to delete contents, "
                    "or find a new dst_dir"
                )

    if not dst_dir.exists():
        dst_dir.mkdir(parents=True)

    host = ftplib.FTP('ftp.ensembl.org')
    host.login()
    release_dir = f'pub/release-{default_ensembl_version}/gff3'
    directory_list = sorted(list(host.nlst(release_dir)))

    if specific_files is not None:
        entry_list = copy.deepcopy(specific_files)
    else:
        entry_list = []

    for ii, ftp_dir in enumerate(directory_list):
        if n_limit is not None and len(entry_list) > n_limit:
            break
        if ii > 0 and ii % 25 == 0:
            print(f"considered {ii} dir; {len(entry_list)} entries")
        file_path_list = list(host.nlst(ftp_dir))
        chosen = None
        for file_path in file_path_list:
            if file_path.endswith(f'{default_ensembl_version}.gff3.gz'):
                assert chosen is None
                chosen = file_path
        if chosen is not None:
            this = {
                'url': f'https://ftp.ensembl.org/{chosen}',
                'assembly_id': 'placeholder.000'
            }
            entry_list.append(this)

    if n_limit is not None:
        print(f"limiting entries to first {n_limit}")
        entry_list = entry_list[:n_limit]
    print(f'processing {len(entry_list)} entries')

    failed_file_lookup = dict()
    serialized_species = set()
    valid_files = []
    n_entries = len(entry_list)
    update_every = max(1, n_entries//10)
    for ii, entry in enumerate(entry_list):
        this_tmp_dir = pathlib.Path(
            tempfile.mkdtemp(dir=tmp_dir)
        )
        if ii > 0 and ii % update_every == 0:
            dur = time.time()-t0
            print(f'processed {ii} of {n_entries} in {dur:.2e} seconds')
        try:
            result_path = serialize_bkbit_gff3(
                content_url=entry['url'],
                assembly_id=entry['assembly_id'],
                dst_dir=dst_dir,
                serialized_species=serialized_species,
                tmp_dir=this_tmp_dir
            )
            valid_files.append(
                {"type": "bkbit",
                 "path": str(result_path.resolve().absolute())
                 }
            )
        except Exception:
            msg = traceback.format_exc()
            failed_file_lookup[entry['url']] = msg
        finally:
            file_utils.clean_up(this_tmp_dir)

    dur = (time.time()-t0)/60.0
    print(f'SUCCESS\nthat took {dur:.2e} minutes')
    with open(failure_log_path, "w") as dst:
        dst.write(json.dumps(failed_file_lookup))

    return valid_files


def serialize_bkbit_gff3(
        content_url,
        assembly_id,
        dst_dir,
        serialized_species,
        tmp_dir=None):

    gff3 = gene_translator.Gff3.from_url(
        content_url=content_url,
        assembly_accession=assembly_id,
        assembly_strain=None,
        log_level="WARNING",
        log_to_file=False,
        use_tqdm=False,
        tmp_dir=tmp_dir
    )
    gff3.parse()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        result = gff3.serialize_to_jsonld()

    taxon = None
    for element in result['@graph']:
        if element['category'][0] == 'biolink:OrganismTaxon':
            taxon = element
    if taxon is None:
        raise RuntimeError(
            f"Could not find taxon for {content_url}"
        )

    taxon_id = taxon['id']
    if taxon_id in serialized_species:
        print(
            f"    == skipping {content_url.split('/')[-1]}; "
            "already serialized"
        )
        raise RuntimeError(f"{taxon_id} already serialized")

    species_name = taxon['full_name'].lower().replace(' ', '_')
    dst_path = dst_dir / f'{species_name}.jsonld'
    if dst_path.exists():
        raise RuntimeError(f"duplicate {dst_path}")

    with open(dst_path, 'w') as dst:
        dst.write(json.dumps(result, indent=2))
    serialized_species.add(taxon_id)
    dst_path = pathlib.Path(dst_path)
    return dst_path
