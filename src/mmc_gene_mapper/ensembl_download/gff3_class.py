"""
Define a subclass of bkbit's Gff3 class with methods patche to fit the needs
of this project.
"""
import pathlib
import re
from urllib.parse import urlparse

from bkbit.models import genome_annotation as genome_annotation
import bkbit.data_translators.genome_annotation_translator as gene_translator


class Gff3Patch(gene_translator.Gff3):
    """
    Subclass of Gff3 that can write jsonld to a dict, rather
    than print it to stdout
    """

    def serialize_to_jsonld(
        self, exclude_none: bool = True, exclude_unset: bool = False
    ):
        """
        Serialize the object and either write it
        to the specified output file or print it to the CLI.

        Parameters:
            exclude_none (bool):
                Whether to exclude None values in the output.
            exclude_unset (bool):
                Whether to exclude unset values in the output.

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
            data.append(
                ck.dict(
                    exclude_none=exclude_none,
                    exclude_unset=exclude_unset
                )
            )
        for ga in self.gene_annotations.values():
            data.append(
                ga.dict(
                    exclude_none=exclude_none,
                    exclude_unset=exclude_unset
                )
            )

        output_data = {
            "@context": (
                "https://raw.githubusercontent.com/brain-bican/"
                "models/main/jsonld-context-autogen/"
                "genome_annotation.context.jsonld"
            ),
            "@graph": data,
        }
        return output_data

    def parse_url(self):
        """
        Parses the content URL and extracts information
        about the genome annotation.

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
        ncbi_pattern = (
            r"/genomes/all/annotation_releases/(\d+)"
            "(?:/(\d+))?/(GCF_\d+\.\d+)[_-]([^/]+)/"  # noqa W605
            "(GCF_\d+\.\d+)[_-]([^/]+)_genomic\.gff\.gz"  # noqa W605
        )

        # Parse the URL to get the path
        parsed_url = urlparse(self.content_url)
        path = parsed_url.path

        # Determine if the URL is from NCBI or Ensembl and extract information
        if "ncbi" in parsed_url.netloc:
            ncbi_match = re.search(ncbi_pattern, path)
            if ncbi_match:
                return {
                    "authority": genome_annotation.AuthorityType.NCBI,
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
                "authority": genome_annotation.AuthorityType.ENSEMBL,
                "release_version": release_version,
                "scientific_name": scientific_name,
                "assembly_name": assembly_name,
            }

        # If no match is found, return None
        return None
