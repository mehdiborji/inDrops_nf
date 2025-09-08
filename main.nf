#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { PREPROCESS_FASTQ } from './modules/fastq_processing'
include { FASTP } from './modules/fastp_qc'
include { STAR_INDEX } from './modules/genome_indexing'
include { GET_INDROPS_V3_WHITELIST } from './modules/whitelist_generation'
include { STARSOLO_ALIGN } from './modules/star_alignment'


def helpMessage() {
    log.info"""
    Usage:

    nextflow run main.nf --input_csv samplesheet.csv --genome_fasta genome.fa --gtf_file genes.gtf -profile full
    
    nextflow run main.nf -profile test

    Required arguments:
    --input_csv         CSV file with sample ID and input FASTQ files  (ID,R1,R2,R3,R4)
    --genome_fasta      Path to genome FASTA file
    --gtf_file          Path to GTF annotation file

    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

if (!params.input_csv) {
    error "Please provide --input_csv parameter"
}
if (!params.genome_fasta) {
    error "Please provide --genome_fasta parameter"
}
if (!params.gtf_file) {
    error "Please provide --gtf_file parameter"
}

workflow {
    // Create input channel from CSV
    // Expected format: sample_id,read1,read2,read3
    input_ch = Channel
        .fromPath(params.input_csv)
        .splitCsv(header: true)
        .map { row -> 
            def sample_id = row.ID
            def read1 = file(row.R1)
            def read2 = file(row.R2)
            def read4 = file(row.R4)
            return [sample_id, read1, read2, read4]
        }
    
    // Process 1: Python preprocessing of 3 input files into one paired-end set
    PREPROCESS_FASTQ(input_ch, params.read_limit)

    // Process 2: fastp quality control and trimming
    FASTP(PREPROCESS_FASTQ.out.reads)

    // Process 3: STAR index building
    STAR_INDEX(file(params.genome_fasta), file(params.gtf_file), params.indexN)

    // Process 4: Making 384*384 16bp whitelist from 384 8bp plate barcodes
    GET_INDROPS_V3_WHITELIST(file(params.whitelist_half))

    // Process 5: STARsolo alignment and quantification
    STARSOLO_ALIGN(FASTP.out.reads,file(params.star_index_dir),GET_INDROPS_V3_WHITELIST.out.whitelist)

}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}