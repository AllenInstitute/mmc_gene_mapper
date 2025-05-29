import json
import mmc_gene_mapper.ensembl_download.scraper as scraper


def main():

    dst_dir = 'output/'

    human = {
        "url": "https://ftp.ensembl.org/pub/release-101/gff3/homo_sapiens/Homo_sapiens.GRCh38.101.gff3.gz",
        "assembly_id": "NCBIAssembly:GCF_000001405.39"
    }

    mouse = {
        "url": "https://ftp.ensembl.org/pub/release-98/gff3/mus_musculus/Mus_musculus.GRCm38.98.gff3.gz",
        "assembly_id": "NCBIAssembly:GCF_000001635.26"
    }

    valid_files = scraper.scrape_ensembl(
        default_ensembl_version=114,
        dst_dir=dst_dir,
        failure_log_path="failures.20250528.json",
        specific_files=[mouse, human]
    )

    with open('valid_files.json', 'w') as dst:
        dst.write(json.dumps(valid_files, indent=2))


if __name__ == "__main__":
    main()
