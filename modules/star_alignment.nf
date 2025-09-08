
process STARSOLO_ALIGN {
    
    tag "$sample_id"

    publishDir params.results_folder, mode: 'symlink'
    
    input:
    tuple val(sample_id), path(read1), path(read2)
    path star_index_dir
    path whitelist
    
    output:
    tuple val(sample_id), path("${sample_id}/"), emit: results_folder

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
        --outSAMtype None \
        --soloOutFormatFeaturesGeneField3 - \
        --soloCellReadStats Standard
    """
    
}
