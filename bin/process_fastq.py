#!/usr/bin/env python
"""
Preprocess reads script - converts 3 input files to 2 paired-end reads
"""

import argparse
import pysam
import logging
import json
from collections import defaultdict

def setup_logging(log_file):
    """Setup logging configuration"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def process_fastq(sample_id, input_r1_fastq, input_r2_fastq, input_r4_fastq, 
                  output_r1_fastq, output_r2_fastq, bc_raw_count_json, limit=None):

    print_sample_reads = False
    
    if limit < 200:
        print_sample_reads = True
    
    R1_clean = open(output_r1_fastq, "w")
    R2_clean = open(output_r2_fastq, "w")
    
    total_reads = 0
    total_failed_reads = 0
    total_accepted_reads = 0
    bc1_bc2_raw_count = defaultdict(int)
    
    with pysam.FastxFile(input_r1_fastq) as R1, \
         pysam.FastxFile(input_r2_fastq) as R2, \
         pysam.FastxFile(input_r4_fastq) as R4:
        
        for r1, r2, r4 in zip(R1, R2, R4):
            total_reads += 1
    
            accept_read = True
    
            rname = r1.name
            
            seq2 = r2.sequence
            qual2 = r2.quality
            
            seq4 = r4.sequence
            qual4 = r4.quality
    
            if print_sample_reads:
                logging.info(f"{seq2} {seq4}")
            
            if len(r1.sequence)<30 or len(seq2)<8 or len(seq4)<14:
                accept_read = False
            
            if seq2.count('G')>6 or seq4.count('G')>12:
                accept_read = False
            
            new_seq2 = seq2 + seq4
            new_qual2 = qual2 + qual4
            
            if accept_read:
    
                total_accepted_reads += 1
    
                bc1_bc2_raw_count[seq2 + seq4[:8]] += 1
    
                R1_clean.write(f"@{rname}\n{r1.sequence}\n+\n{r1.quality}\n")
                R2_clean.write(f"@{rname}\n{new_seq2}\n+\n{new_qual2}\n")
                
                if print_sample_reads:
                    logging.info(f"{new_seq2} \n")
                    
            else:
                total_failed_reads += 1
            
                if print_sample_reads:
                    logging.info("Failed reads \n")
            
            if total_reads % 5000000 == 0:
                logging.info(f"Processed {total_reads} reads")
            if limit and total_reads > limit:
                break
    R1_clean.close()
    R2_clean.close()

    logging.info(f"Processed total {total_reads} reads")
    logging.info(f"{total_failed_reads} failed reads")
    logging.info(f"{total_accepted_reads} accepted reads")
    
    with open(bc_raw_count_json, "w") as json_file:
        json.dump(bc1_bc2_raw_count, json_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Preprocess 3 input files to 2 paired-end reads')
    parser.add_argument('--sample_id', required=True, help='Sample identifier')
    parser.add_argument('--input_r1_fastq', required=True, help='Input Read1 FASTQ')
    parser.add_argument('--input_r2_fastq', required=True, help='Input Read2 FASTQ')
    parser.add_argument('--input_r4_fastq', required=True, help='Input Read4 FASTQ')
    parser.add_argument('--output_r1_fastq', required=True, help='Output R1 FASTQ')
    parser.add_argument('--output_r2_fastq', required=True, help='Output R2 FASTQ')
    parser.add_argument('--bc_raw_count_json', required=True, help='Raw Barcode Count JSON')
    parser.add_argument('--limit', type=int, default=1e9, help='Number of reads to test')
    parser.add_argument('--log', required=True, help='Log file')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log)
    
    # Run preprocessing
    process_fastq(
        args.sample_id,
        args.input_r1_fastq,
        args.input_r2_fastq,
        args.input_r4_fastq,
        args.output_r1_fastq,
        args.output_r2_fastq,
        args.bc_raw_count_json,
        args.limit
    )
    
    logging.info("Preprocessing completed successfully")