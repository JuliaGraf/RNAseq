# RNAseq

## Usage

First, prepare a sample datasheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,strandedness
SHAM_OXY,path_to_fastq_1.fastq.gz,path_to_fastq_2.fastq.gz,auto
SNI_OXY,path_to_fastq_1.fastq.gz,path_to_fastq_2.fastq.gz,auto
```

Now, you can run the pipeline using:

```bash
nextflow run ./workflows/main.nf --input <path/to/input_samplesheet>  --genomeFasta <path/to/input_genomeFasta> --gtfFile <path/to/gtfFile>  -entry RNASEQ
```

> [!NOTE]
> To run the pipeline you need to download the genome from [[NCBI]](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz) and unzip the gtf and genome file. 
> The input read files also have to be unzipped.
