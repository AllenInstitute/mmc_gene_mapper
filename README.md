# MapMyCells Gene Mapper

## Overview

This code provides a tool to map genes from one species/authority (where
by "authority" we mean an institution like ENSEMBL or NCBI) to another
species/authority using documented cross-authority and crosss-species
(orthologous) gene equvialencies. This code was originally written to allow
users to map data collected from an arbitrary species to the species-specific
cell type taxonomies supported by the Allen Institute's
[MapMyCells](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells)
tool and
[cell_type_mapper library](https://github.com/AllenInstitute/cell_type_mapper).

It can, in principle, be applied to other use cases where it would be useful
to map genes from one identification scheme to another.

This is a re-implementation in Python of the functionality available through the R [GeneOrthology](https://github.com/AllenInstitute/GeneOrthology) package.

## Installation

To install this library, either

A) Clone the repository and run `pip install .` from the root directory of the
repository

or

B) Run
```
pip install "mmc_gene_mapper @ git+https://github.com/AllenInstitute/mmc_gene_mapper"
```

To install a specific version, you can run
```
pip install "mmc_gene_mapper @ git+https://github.com/AllenInstitute/mmc_gene_mapper@{VERSION}"
```
where `{VERSION}` is one of the
[valid tags listed in this repository.](https://github.com/AllenInstitute/mmc_gene_mapper/tags)

## Use

[This Jupyter notebook](https://github.com/AllenInstitute/mmc_gene_mapper/blob/main/notebooks/gene_mapper_demo.ipynb)
demonstrates how to use the code in this repository.

In broad strokes, you must first create a sqlite database file containing the valid gene mappings as determined from data published by NCBI and ENSEMBL. This can be done either programmatically, according to cell [5] the notebook
referenced above, or using the command line tool
```
python -m mmc_gene_mapper.cli.create_db_file --help
```
The `--help` will cause the arguments expected by the command line tool
to appear in stdout. A good default would be
```
python -m mmc_gene_maper.cli.create_db_file \
--db_path data/my_gene_mapper.db \
--local_dir data/downloaded_files/ \
--ensembl_version 114
```

**Note:** This will produce a 15 GB file on disk.

Once that file has been created, you can instantiate the class
`MMCGeneMapper` (again, see cell [5] of the exampe notebook), and you are
ready to map genes.

## Level of support

We are providing this tool to the community and any and all who want to use it.
Issues and pull requests are welcome, however, this code is also intended
as part of the backend for the Allen Institute Brain Knowledge Platform. As
such, issues and pull requests may be declined if they interfere with
the functionality required to support that service.
