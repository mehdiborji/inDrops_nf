#!/usr/bin/env python
"""
Process gel_barcode2_list to get the 384*384 possible barcodes
"""

import argparse
from itertools import product
import pandas as pd

def rev(seq):
    trans = bytes.maketrans(b"ATGC", b"TACG")
    return seq.encode().translate(trans)[::-1].decode()

def process_plate_bcs(gel_barcode2_list, output_whitelist):

    gel_bc2 = pd.read_csv(gel_barcode2_list, header=None, names=["seq"])
    
    gel_bc2['rev_seq'] = gel_bc2.map(rev)
    
    v3_product = product(gel_bc2['seq'],gel_bc2['rev_seq'])
    
    v3_barcodes = [bc + rev for bc, rev in v3_product]

    pd.Series(v3_barcodes).to_csv(output_whitelist, header=None, index=None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process gel_barcode2_list')
    parser.add_argument('--gel_barcode2_list', required=True, help='384 well plate 8bp list')
    parser.add_argument('--output_whitelist', required=True, help='final combinations whitelist')

    args = parser.parse_args()

    process_plate_bcs(
        args.gel_barcode2_list,
        args.output_whitelist
    )
