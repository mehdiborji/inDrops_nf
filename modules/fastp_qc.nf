
/*
QC Filtering and Adapter Trimming with FASTP
Length 22 ensures R2 has 16nt for BC and 6nt for UMI
*/

process FASTP {
    
    tag "$sample_id"
    publishDir params.results_folder, mode: 'symlink'
    
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