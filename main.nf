#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
//include { PREPROCESS } from './modules/preprocess'
//include { FASTP } from './modules/fastp'
//include { STARSOLO } from './modules/starsolo'

// Help message
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --input_csv samplesheet.csv --genome_fasta genome.fa --gtf_file genes.gtf [options]
    
    Required arguments:
    --input_csv         CSV file with sample information (sample_id,file1,file2,file3)
    --genome_fasta      Path to genome FASTA file
    --gtf_file          Path to GTF annotation file
    
    Optional arguments:
    --star_index        Path to pre-built STAR index (if not provided, will be built)
    --outdir            Output directory (default: results)
    
    Example:
    nextflow run main.nf \\
        --input_csv samplesheet.csv \\
        --genome_fasta genome.fa \\
        --gtf_file genes.gtf \\
        --outdir results
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Validate required parameters
if (!params.input_csv) {
    error "Please provide --input_csv parameter"
}


process PREPROCESS {
    
    tag "$sample_id"
    publishDir params.results_folder, mode: 'symlink'
    
    input:
    tuple val(sample_id), path(read1), path(read2), path(read4)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq"), emit: reads
    path "${sample_id}_preprocessing.log", emit: log
    
    script:
    """
    preprocess.py \\
        --sample_id ${sample_id} \\
        --input_r1_fastq ${read1} \\
        --input_r2_fastq ${read2} \\
        --input_r4_fastq ${read4} \\
        --output_r1_fastq ${sample_id}_R1.fastq \\
        --output_r2_fastq ${sample_id}_R2.fastq \\
        --bc_raw_count_json ${sample_id}_raw_bcs.json \\
        --limit 2000000 \\
        --log ${sample_id}_preprocessing.log
    """

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
    
    // Process 1: Python preprocessing
    // Takes tuple of 3 files, outputs 2 paired-end reads
    PREPROCESS(input_ch)
}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}