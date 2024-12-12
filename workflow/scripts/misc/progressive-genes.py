#!/usr/bin/env python3

import argparse
import pandas as pd
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description='Find genes/isoforms that are differentially expressed in the same direction across comparisons'
    )
    parser.add_argument('--comps', nargs='+', required=True,
                       help='Comparisons in format control_intermediate_resistant')
    parser.add_argument('--padj-threshold', type=float, default=0.05,
                       help='Adjusted p-value threshold')
    parser.add_argument('--fc-threshold', type=float, default=2.0,
                       help='Fold change threshold')
    parser.add_argument('--gene-level', action='store_true',
                       help='Analyze gene-level differential expression')
    parser.add_argument('--isoform-level', action='store_true',
                       help='Analyze isoform-level differential expression')
    return parser.parse_args()

def process_comparison(sus, intermediate, res, padj_threshold, upper_fc, lower_fc, level='gene'):
    """Process a single comparison for either genes or isoforms."""
    
    # Set paths and column names based on analysis level
    if level == 'gene':
        base_path = 'results/genediff'
        id_col = 'GeneID'
        padj_col = 'padj'
    else:  # isoform
        base_path = 'results/isoformdiff'
        id_col = 'TranscriptID'
        padj_col = 'qval'
    
    # Read the comparison files
    one = pd.read_csv(f"{base_path}/{intermediate}_{res}.csv")
    two = pd.read_csv(f"{base_path}/{sus}_{intermediate}.csv")
    
    # Filter for up and down regulated genes/isoforms
    up1 = one[
        (one['FC'] > upper_fc) & 
        (one[padj_col] < padj_threshold)
    ]
    down1 = one[
        (one['FC'] < lower_fc) & 
        (one[padj_col] < padj_threshold)
    ]
    
    up2 = two[
        (two['FC'] > upper_fc) & 
        (two[padj_col] < padj_threshold)
    ]
    down2 = two[
        (two['FC'] < lower_fc) & 
        (two[padj_col] < padj_threshold)
    ]
    
    # Find intersections
    intersect_up = pd.merge(
        up1, up2,
        on=id_col,
        suffixes=('_field', '_lab')
    )
    intersect_down = pd.merge(
        down1, down2,
        on=id_col,
        suffixes=('_field', '_lab')
    )
    
    # Save results
    output_base = f"{base_path}/{sus}_{intermediate}_{res}"
    intersect_up.to_csv(f"{output_base}.up.progressive.tsv", sep='\t', index=False)
    intersect_down.to_csv(f"{output_base}.down.progressive.tsv", sep='\t', index=False)

def main():
    args = parse_args()
    
    # Calculate lower fold change threshold
    upper_fc = args.fc_threshold
    lower_fc = 1.0 / upper_fc
    
    # Process each comparison
    for comp in args.comps:
        sus, intermediate, res = comp.split('_')
        
        if args.gene_level:
            process_comparison(
                sus, intermediate, res,
                args.padj_threshold, upper_fc, lower_fc,
                level='gene'
            )
            
        if args.isoform_level:
            process_comparison(
                sus, intermediate, res,
                args.padj_threshold, upper_fc, lower_fc,
                level='isoform'
            )

if __name__ == '__main__':
    main()