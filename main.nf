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
if (!params.genome_fasta) {
    error "Please provide --genome_fasta parameter"
}
if (!params.gtf_file) {
    error "Please provide --gtf_file parameter"
}

process PREPROCESS_FASTQ {
    
    tag "$sample_id"
    publishDir params.results_folder, mode: 'move'
    
    input:
    tuple val(sample_id), path(read1), path(read2), path(read4)
    val(limit)
    
    output:
    tuple val(sample_id), path("${sample_id}_R1.fastq"), path("${sample_id}_R2.fastq"), emit: reads
    path "${sample_id}_preprocessing.log", emit: log
    
    script:
    """
    process_fastq.py \
        --sample_id ${sample_id} \
        --input_r1_fastq ${read1} \
        --input_r2_fastq ${read2} \
        --input_r4_fastq ${read4} \
        --output_r1_fastq ${sample_id}_R1.fastq \
        --output_r2_fastq ${sample_id}_R2.fastq \
        --bc_raw_count_json ${sample_id}_raw_bcs.json \
        --limit ${limit} \
        --log ${sample_id}_preprocessing.log
    """

}

process FASTP {
    
    tag "$sample_id"
    publishDir params.results_folder, mode: 'move'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    
    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq"), path("${sample_id}_trimmed_R2.fastq"), emit: reads
    path "${sample_id}_fastp.json", emit: json
    path "${sample_id}_fastp.html", emit: html
    
    script:
    """
    fastp \
        --in1 ${read1} \
        --in2 ${read2} \\
        --out1 ${sample_id}_trimmed_R1.fastq \
        --out2 ${sample_id}_trimmed_R2.fastq \\
        --json ${sample_id}_fastp.json \
        --html ${sample_id}_fastp.html \\
        --thread ${task.cpus} \
        --length_required 22 \\
        --detect_adapter_for_pe
    """
}


process STAR_INDEX {
    
    container 'longread_nf'

    tag "STAR_INDEX"
    publishDir params.star_index_dir, mode: 'move'
    
    input:
    path genome_fasta
    path gtf_file
    
    output:
    path ".", emit: index_dir
    
    when:
    // Only run if index directory doesn't exist or is empty
    def indexDir = file(params.star_index_dir)
    !indexDir.exists() || indexDir.list().size() == 0

    script:
    """
    STAR \
        --runMode genomeGenerate \\
        --genomeDir . \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${gtf_file} \\
        --runThreadN ${task.cpus} \\
        --sjdbOverhang 60 \\
        --genomeSAindexNbases 14
    """

}

process GET_WHITELIST {
    
    publishDir params.results_folder, mode: 'move'
    
    input:
    path bc_list
    
    output:
    path("output_whitelist"), emit: whitelist
    
    script:
    """
    process_plate_bcs.py \\
        --gel_barcode2_list ${bc_list} \\
        --output_whitelist output_whitelist
    """

}


process STARSOLO_ALIGN {
    
    container 'longread_nf'

    tag "STARSOLO_ALIGN_$sample_id"

    publishDir params.results_folder, mode: 'move'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    path star_index_dir
    path whitelist
    
    output:
    tuple val(sample_id), path("${sample_id}/"), emit: results
    tuple val(sample_id), path("${sample_id}/*Aligned.sortedByCoord.out.bam"), emit: bam
    
    script:

    """
    STAR \
        --runThreadN ${task.cpus} \
        --readFilesIn  ${read1} ${read2} \
        --genomeDir ${star_index_dir} \
        --outFileNamePrefix ${sample_id}/ \
        --soloCBwhitelist ${whitelist} \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 16 \
        --soloUMIstart 17 \
        --soloUMIlen 6 \
        --soloFeatures Gene Velocyto GeneFull \
        --soloUMIdedup Exact \
        --soloMultiMappers EM \
        --soloCellFilter None \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes AS CR UR CB UB GX GN \
        --soloOutFormatFeaturesGeneField3 - \
        --soloCellReadStats Standard
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
    PREPROCESS_FASTQ(input_ch, 2000000000)

    //PREPROCESS.out.reads.view()
    // Process 2: fastp quality control and trimming
    // Takes paired-end reads, outputs trimmed reads
    FASTP(PREPROCESS_FASTQ.out.reads)

    // Process 3: STAR index building (conditional)
    // Check if index exists, build if not available
    STAR_INDEX(file(params.genome_fasta), file(params.gtf_file))

    // Process 4: make 384*384 16bp whitelist from 384 8bp plate barcodes
    GET_WHITELIST(file(params.whitelist_half))

    // Process 5: STARsolo alignment and quantification
    // Takes trimmed reads and STAR index, outputs alignment and count matrix
    STARSOLO_ALIGN(FASTP.out.reads,file(params.star_index_dir),GET_WHITELIST.out.whitelist)

}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}