/*
Generate genome index from fasta and gtf file
storeDir directive will serve as permenant cache
Depending on genome being small (Mito genome) or whole genome
genomeSAindexNbases needs to be dynamically modified
*/

process STAR_INDEX {
    
    tag "STAR_INDEX"
    storeDir params.star_index_dir
    
    input:
    path genome_fasta
    path gtf_file
    val indexN
    
    output:
    path "*", emit: index_dir
    
    script:
    """
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir . \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile ${gtf_file} \\
        --runThreadN ${task.cpus} \\
        --sjdbOverhang 60 \\
        --genomeSAindexNbases ${indexN}
    """

}