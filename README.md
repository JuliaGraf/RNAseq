# RNAseq

## Usage

First, prepare a sample datasheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness
sample1,<path/to/fastq_1_file>,<path/to/fastq_2_file>
sample2,<path/to/fastq_1_file>,<path/to/fastq_2_file>
```

Now, you can run the pipeline using:

```bash
nextflow run ./workflows/main.nf \
  --input <path/to/input_samplesheet> \
  --genomeFasta <path/to/input_genomeFasta> \
  --gtfFile <path/to/gtfFile> \
  [--amd/--arm] \
  -entry RNASEQ
```
> [!IMPORTANT]
> The pipeline can be run with either `--amd` or `--arm`, depending on the architecture of the user's operating system. If nothing is specified, the pipeline is run using `--amd` by default.

> [!NOTE]
> To run the pipeline on the test data, you need to download the [gtf](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz) and [genome](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz) files and unzip them. 
