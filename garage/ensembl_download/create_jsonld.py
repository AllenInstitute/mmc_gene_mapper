import json
import pandas as pd
import pathlib
import time
import warnings

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
        Serialize the object and either write it to the specified output file or print it to the CLI.

        Parameters:
            exclude_none (bool): Whether to exclude None values in the output.
            exclude_unset (bool): Whether to exclude unset values in the output.

        Returns:
            None
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


def main():
    df = pd.read_csv(
        'data/20240125_gars_genome_annotation_tracking.csv'
    )

    dst_dir = pathlib.Path('output')
    assert dst_dir.is_dir()

    t0 = time.time()
    df = df[df.genome_annotation_content_url.str.contains('ensembl.org')]
    this_pass = set()
    for entry in df.iterrows():
        entry = entry[1]
        species = entry.organism_taxon_full_name.lower()
        print(f'PROCESSING {species} {time.time()-t0:.2e}')
        species = species.replace(' ', '_')
        dst_path = dst_dir / f'{species}.jsonld'
        if dst_path in this_pass:
            print(f'    duplicate {dst_path}')
        this_pass.add(dst_path)
        parse_gff3(entry, dst_path=dst_path)

    dur = (time.time()-t0)/60.0
    print(f'that took {dur:.2e} minutes')

def parse_gff3(entry, dst_path):
    content_url = entry.genome_annotation_content_url
    assembly_id = entry.genome_assembly_id

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

    with open(dst_path, 'w') as dst:
        dst.write(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
