import ftplib
import json
import pandas as pd
import pathlib
import time
from urllib.parse import urlparse
import warnings

from bkbit.models import genome_annotation as ga
import bkbit.data_translators.genome_annotation_translator as gene_translator

import mmc_gene_mapper.utils.file_utils as file_utils


class Gff3Patch(gene_translator.Gff3):
    """
    Subclass of Gff3 that can write jsonld to a dict, rather
    than print it to stdout
    """

    def serialize_to_jsonld(
        self, exclude_none: bool = True, exclude_unset: bool = False
    ):
        """
        Serialize the object and either write it to the specified output file or print it to the CLI.

        Parameters:
            exclude_none (bool): Whether to exclude None values in the output.
            exclude_unset (bool): Whether to exclude unset values in the output.

        Returns:
            Dict

        Notes
        -----
        Overrode default implementation so that we can return the dict
        containing the jsonld graph
        """

        data = [
            self.organism_taxon.dict(
                exclude_none=exclude_none, exclude_unset=exclude_unset
            ),
            self.genome_assembly.dict(
                exclude_none=exclude_none, exclude_unset=exclude_unset
            ),
            self.genome_annotation.dict(
                exclude_none=exclude_none, exclude_unset=exclude_unset
            ),
        ]
        for ck in self.checksums:
            data.append(ck.dict(exclude_none=exclude_none, exclude_unset=exclude_unset))
        for ga in self.gene_annotations.values():
            data.append(ga.dict(exclude_none=exclude_none, exclude_unset=exclude_unset))

        output_data = {
            "@context": "https://raw.githubusercontent.com/brain-bican/models/main/jsonld-context-autogen/genome_annotation.context.jsonld",
            "@graph": data,
        }
        return output_data


    def parse_url(self):
        """
        Parses the content URL and extracts information about the genome annotation.

        Returns:
            A dictionary containing the following information:
            - 'authority': The authority type (NCBI or ENSEMBL).
            - 'taxonid': The taxon ID of the genome.
            - 'release_version': The release version of the genome annotation.
            - 'assembly_accession': The assembly accession of the genome.
            - 'assembly_name': The name of the assembly.
            - 'species': The species name (only for ENSEMBL URLs).

        Notes
        -----
        Overrode default implementation so that we can parse URLs with
        non-default file name schemas
        """
        # Define regex patterns for NCBI and Ensembl URLs
        # NCBI : [assembly accession.version]_[assembly name]_[content type].[optional format]
        # ENSEMBL :  <species>.<assembly>.<_version>.gff3.gz -> organism full name, assembly name, genome version
        ncbi_pattern = r"/genomes/all/annotation_releases/(\d+)(?:/(\d+))?/(GCF_\d+\.\d+)[_-]([^/]+)/(GCF_\d+\.\d+)[_-]([^/]+)_genomic\.gff\.gz"
        ensembl_pattern = (
            r"/pub/release-(\d+)/gff3/([^/]+)/([^/.]+)\.([^/.]+)\.([^/.]+)\.gff3\.gz"
        )

        # Parse the URL to get the path
        parsed_url = urlparse(self.content_url)
        path = parsed_url.path

        # Determine if the URL is from NCBI or Ensembl and extract information
        if "ncbi" in parsed_url.netloc:
            ncbi_match = re.search(ncbi_pattern, path)
            if ncbi_match:
                return {
                    "authority": ga.AuthorityType.NCBI,
                    "taxonid": ncbi_match.group(1),
                    "release_version": (
                        ncbi_match.group(2)
                        if ncbi_match.group(2)
                        else ncbi_match.group(4)
                    ),
                    "assembly_accession": ncbi_match.group(3),
                    "assembly_name": ncbi_match.group(6),
                }

        elif "ensembl" in parsed_url.netloc:
            file_path = pathlib.Path(parsed_url.path).name
            file_path = file_path.replace('.gff3.gz', '')

            scientific_name = file_path.split('.')[0]
            file_path = file_path.replace(scientific_name+'.', '')
            release_version = file_path.split('.')[-1]
            file_path = file_path.replace('.' + release_version, '')
            assembly_name = file_path

            return {
                "authority": ga.AuthorityType.ENSEMBL,
                "release_version": release_version,
                "scientific_name": scientific_name,
                "assembly_name":assembly_name,
            }

        # If no match is found, return None
        return None



def main():
    #df = pd.read_csv(
    #    'data/20240125_gars_genome_annotation_tracking.csv'
    #)

    dst_dir = pathlib.Path('output')

    t0 = time.time()
    host = ftplib.FTP('ftp.ensembl.org')
    host.login()
    release_dir = 'pub/release-114/gff3'
    directory_list = sorted(list(host.nlst(release_dir)))
    failed_file_list = []
    for entry in directory_list:
        print(entry,time.time()-t0)
        file_path_list = list(host.nlst(entry))
        chosen = None
        for file_path in file_path_list:
            if file_path.endswith('114.gff3.gz'):
                assert chosen is None
                chosen = file_path
        assert chosen is not None
        try:
            serialize_bkbit_gff3(
                content_url=f'https://ftp.ensembl.org/{chosen}',
                assembly_id='placeholder.000',
                dst_dir=dst_dir
            )
        except Exception:
            print(f'    {chosen} failed')
            failed_file_list.append(chosen)
    print(f'failed files\n{json.dumps(failed_file_list, indent=2)}')
    dur = (time.time()-t0)/60.0
    print(f'SUCCESS\nthat took {dur:.2e} minutes')

def serialize_bkbit_gff3(
        content_url,
        assembly_id,
        dst_dir):

    gff3 = Gff3Patch(
        content_url=content_url,
        assembly_accession=assembly_id,
        assembly_strain=None,
        log_level="WARNING",
        log_to_file=False
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
    species_name = taxon['full_name'].lower().replace(' ', '_')
    dst_path = dst_dir / f'{species_name}.jsonld'
    if dst_path.exists():
        print(f"duplicate {dst_path}")

    with open(dst_path, 'w') as dst:
        dst.write(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
