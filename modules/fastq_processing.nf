/*
Take 3 or 4 input fastq files with inDropsV3 format
Discard short and polyG reads
Concatenate R2 and R4 which contain Cell Barcode and UMI
*/

process PREPROCESS_FASTQ {
    
    tag "$sample_id"
    publishDir params.results_folder, mode: 'symlink'
    
    input:
    tuple val(sample_id), path(read1), path(read2), path(read4)
    val read_limit
    
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
        --limit ${read_limit} \
        --log ${sample_id}_preprocessing.log
    """

}