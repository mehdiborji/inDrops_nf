# inDrops v3 Bioinformatics Pipeline and Downstream Analysis

This pipeline only leverages a tacit knowledge about the library read structure of inDrops v3 assay provided in the function [`stable_barcode_names`](https://github.com/indrops/indrops/blob/master/indrops.py#L352-L381), where the inDrops v3 barcodes are defined.


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
│   ├── process_fastq.py
│   └── process_plate_bcs.py
├── conda.yml
├── data
│   ├── gel_barcode2_list.txt
│   ├── gencode_v44_human_chrM
│   │   ├── chrM.fa
│   │   └── chrM.gtf
│   ├── SRA_indrops_sample_sheet_chrM_50k.csv
│   └── SRA_indrops_sample_sheet.csv
├── Dockerfile
├── fastq_chrM_50k
│   ├── SRR16287104_1.fastq.gz
│   ├── SRR16287104_2.fastq.gz
│   ├── SRR16287104_3.fastq.gz
│   ├── SRR16287104_4.fastq.gz
│   ├── SRR16287105_1.fastq.gz
│   ├── SRR16287105_2.fastq.gz
│   ├── SRR16287105_3.fastq.gz
│   ├── SRR16287105_4.fastq.gz
│   ├── SRR16287106_1.fastq.gz
│   ├── SRR16287106_2.fastq.gz
│   ├── SRR16287106_3.fastq.gz
│   ├── SRR16287106_4.fastq.gz
│   ├── SRR16287107_1.fastq.gz
│   ├── SRR16287107_2.fastq.gz
│   ├── SRR16287107_3.fastq.gz
│   └── SRR16287107_4.fastq.gz
├── main.nf
├── modules
│   ├── fastp_qc.nf
│   ├── fastq_processing.nf
│   ├── genome_indexing.nf
│   ├── star_alignment.nf
│   └── whitelist_generation.nf
├── nextflow.config
└── README.md
```
A small dataset of only Mito reads from single-cell experiments each with 4 reads is included. Along with that appropriate sample sheet `csv`, genome `fasta` and annotations `gtf` are included.
Within the `inDrops_nf` directory, build the Docker image:

```
docker built -t indrops_starsolo_nf .
```


## Run the Pipeline

 To submit the sample dataset to Nextflow engine and quantifying the reads, while you are within `inDrops_nf` folder the pipeline can be evoked in `test` mode using the following command.

For binned outputs:

```
nextflow run main.nf -profile test
```

This is what the command output should produce on the screen
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

The outputs are visible in `$HOME/indrops_outputs`.
In addition to intermediate files and reports, for each sample a final folder containing single-cell count information should be available with the following contents:

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