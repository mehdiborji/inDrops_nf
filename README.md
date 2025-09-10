# inDrops v3 Bioinformatics Pipeline and Downstream Analysis

This pipeline only leverages a tacit knowledge about the library read structure of inDrops v3 assay provided in the function [`stable_barcode_names`](https://github.com/indrops/indrops/blob/master/indrops.py#L352-L381), where the inDrops v3 barcodes are defined.

This pipeline is built around [`STARsolo`](https://github.com/alexdobin/STAR) for quantifying inDrops v3 single-cell libraries. The libraries consist of separate read files where `R1` is cDNA, `R2` is the 8 bp barcode, `R3` is the library index, and `R4` is the second barcode plus a 6 bp UMI. To use `STARsolo`, R2 and R4 are concatenated to create a new read. `fastp` is used for QC and adapter trimming.

The input parameters to the pipeline are a genome `FASTA` file and an annotation `GTF` file, as well as a sample sheet in `CSV` format that contains one library per row and four additional columns with paths to R1–R4 for each library. A whitelist of 16 bp barcodes is provided as part of the input and is generated on the fly using the 384-well plate barcodes. Each barcode is corrected within one Hamming distance of this set.

The outputs are count matrices in different formats, including Gene (exonic reads), GeneFull (exonic + intronic reads), and Velocyto (spliced + unspliced reads).


## Install Prerequisites

Docker and Nextflow are required before running the pipeline:
```
https://docs.docker.com/get-started/get-docker/
```
```
https://nextflow.io/docs/latest/install.html
```
## Set Up the Pipeline

Download Nextflow Pipeline
```
cd $HOME
git clone https://github.com/mehdiborji/inDrops_nf
cd inDrops_nf
```

Using `tree` you should be able to see the following:
```
.
├── bin
│   ├── process_fastq.py
│   └── process_plate_bcs.py
├── conda.yml
├── data
│   ├── fastq_chrM_50k
│   │   ├── SRR16287104_1.fastq.gz
│   │   ├── SRR16287104_2.fastq.gz
│   │   ├── SRR16287104_3.fastq.gz
│   │   ├── SRR16287104_4.fastq.gz
│   │   ├── SRR16287105_1.fastq.gz
│   │   ├── SRR16287105_2.fastq.gz
│   │   ├── SRR16287105_3.fastq.gz
│   │   ├── SRR16287105_4.fastq.gz
│   │   ├── SRR16287106_1.fastq.gz
│   │   ├── SRR16287106_2.fastq.gz
│   │   ├── SRR16287106_3.fastq.gz
│   │   ├── SRR16287106_4.fastq.gz
│   │   ├── SRR16287107_1.fastq.gz
│   │   ├── SRR16287107_2.fastq.gz
│   │   ├── SRR16287107_3.fastq.gz
│   │   └── SRR16287107_4.fastq.gz
│   ├── gel_barcode2_list.txt
│   ├── gencode_v44_human_chrM
│   │   ├── chrM.fa
│   │   └── chrM.gtf
│   ├── SRA_indrops_sample_sheet_chrM_50k.csv
│   └── SRA_indrops_sample_sheet.csv
├── Dockerfile
├── main.nf
├── modules
│   ├── fastp_qc.nf
│   ├── fastq_processing.nf
│   ├── genome_indexing.nf
│   ├── star_alignment.nf
│   └── whitelist_generation.nf
├── nextflow.config
└── README.md
```

A small dataset of 50k mitochondrial reads (`data/fastq_chrM_50k`) from a set of libraries in an [inDrops v3 study](https://www.ncbi.nlm.nih.gov/sra?term=SRP340747), each with 4 reads, is included. Along with this, an appropriate sample sheet (`data/SRA_indrops_sample_sheet_chrM_50k.csv`), genome (`data/gencode_v44_human_chrM/chrM.fa`), and annotations (`data/gencode_v44_human_chrM/chrM.fa`) are also provided.


Before running the pipeline, navigate to the `inDrops_nf` directory and build the Docker image using the provided `conda.yml` and `Dockerfile`.
```
docker built -t indrops_starsolo_nf .
```


## Run the Pipeline

To submit the sample dataset to the Nextflow engine and quantify the reads, while in the `inDrops_nf` folder, the pipeline can be invoked in `test` mode using the following command:

```
nextflow run main.nf -profile test
```

The command output should look like this on the screen:

```
 N E X T F L O W   ~  version 25.04.6

Launching `main.nf` [prickly_bassi] DSL2 - revision: 3a9d0a2c71

WARN: Access to undefined parameter `help` -- Initialise it to a default value eg. `params.help = some_value`
executor >  local (14)
[9d/06c1b0] process > PREPROCESS_FASTQ (SRR16287106) [100%] 4 of 4 ✔
[58/bfe8c6] process > FASTP (SRR16287106)            [100%] 4 of 4 ✔
[6d/ed5639] process > STAR_INDEX (STAR_INDEX)        [100%] 1 of 1 ✔
[e2/1d45df] process > GET_INDROPS_V3_WHITELIST       [100%] 1 of 1 ✔
[7c/49be06] process > STARSOLO_ALIGN (SRR16287106)   [100%] 4 of 4 ✔
Pipeline completed at: 2025-09-08T21:07:43.753800-04:00
Execution status: OK
Execution duration: 22.1s
```

The outputs are visible in `$HOME/indrops_outputs`, with symlinks to the original files in the Nextflow work directory, which by default is created in the location where the pipeline was executed.

In addition to intermediate files and reports, a final folder containing single-cell count information should be available for each sample, with the following contents:

```
├── SRR16287105
│   ├── Log.final.out
│   ├── Log.out
│   ├── Log.progress.out
│   ├── SJ.out.tab
│   └── Solo.out
│       ├── Barcodes.stats
│       ├── Gene
│       │   ├── CellReads.stats
│       │   ├── Features.stats
│       │   ├── raw
│       │   │   ├── barcodes.tsv
│       │   │   ├── features.tsv
│       │   │   ├── matrix.mtx
│       │   │   └── UniqueAndMult-EM.mtx
│       │   └── Summary.csv
│       ├── GeneFull
│       │   ├── CellReads.stats
│       │   ├── Features.stats
│       │   ├── raw
│       │   │   ├── barcodes.tsv
│       │   │   ├── features.tsv
│       │   │   ├── matrix.mtx
│       │   │   └── UniqueAndMult-EM.mtx
│       │   └── Summary.csv
│       └── Velocyto
│           ├── Features.stats
│           ├── raw
│           │   ├── ambiguous.mtx
│           │   ├── barcodes.tsv
│           │   ├── features.tsv
│           │   ├── spliced.mtx
│           │   └── unspliced.mtx
│           └── Summary.csv

```