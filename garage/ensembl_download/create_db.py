import pathlib
import mmc_gene_mapper.mapper.mapper as mapper


data_file_spec = []

bkbit_dir = pathlib.Path('output')
bkbit_path_list = [n for n in bkbit_dir.rglob('**/*jsonld')]

data_file_spec = [
    {"type": "bkbit",
     "path": str(pth)}
    for pth in bkbit_path_list
    if 'saccharomyces_cerevisiae' not in pth.name
]

data_file_spec.append(
    {"type": "hmba_orthologs",
     "path": "../../data/db_creation_data/all_gene_ids.csv",
     "name": "HMBA",
     "baseline_species": "human"}
)

print(data_file_spec)

db_path = 'gene_mapper.full_ensembl.db'

gene_mapper = mapper.MMCGeneMapper(
    db_path=db_path,
    local_dir='../../data',
    data_file_spec=data_file_spec,
    clobber=True,
    force_download=False,
    suppress_download_stdout=False
)

print("SUCCESS")
