/*
Generare 384*384 product from 384 plate barcodes
This 144k list of 16bp barcodes will be used as whitelist in STARsolo
*/

process GET_INDROPS_V3_WHITELIST {
    
    publishDir params.results_folder, mode: 'symlink'
    
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