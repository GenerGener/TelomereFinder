#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import re

class TelomereFinder:
    def __init__(self, kmer_sizes, end_region_size=1000):
        """
        Initialize TelomereFinder with k-mer sizes and region size for end analysis
        
        Args:
            kmer_sizes (list): List of k-mer sizes to analyze
            end_region_size (int): Size of region to consider as chromosome end
        """
        self.kmer_sizes = kmer_sizes
        self.end_region_size = end_region_size
        self.kmer_counts = defaultdict(int)
        self.kmer_positions = defaultdict(list)
        
    def process_sequence(self, seq_id, sequence):
        """Process a single sequence for k-mer occurrences"""
        sequence = str(sequence).upper()
        seq_len = len(sequence)
        
        for k in self.kmer_sizes:
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                rev_kmer = str(Seq(kmer).reverse_complement())
                
                # Record position and orientation
                self.kmer_counts[kmer] += 1
                self.kmer_positions[seq_id].append({
                    'kmer': kmer,
                    'start': i,
                    'end': i + k,
                    'orientation': '+',
                    'is_end': i < self.end_region_size or i > seq_len - self.end_region_size
                })
                
                # Also count reverse complement occurrences
                self.kmer_counts[rev_kmer] += 1
    
    def write_bed(self, output_path):
        """Write results in BED format (0-based)"""
        with open(output_path, 'w') as f:
            for seq_id, positions in self.kmer_positions.items():
                for pos in positions:
                    f.write(f"{seq_id}\t{pos['start']}\t{pos['end']}\t"
                           f"{pos['kmer']}\t0\t{pos['orientation']}\n")
    
    def write_gff(self, output_path):
        """Write results in GFF format (1-based)"""
        with open(output_path, 'w') as f:
            f.write("##gff-version 3\n")
            for seq_id, positions in self.kmer_positions.items():
                for pos in positions:
                    f.write(f"{seq_id}\tkmer_finder\ttelomeric_repeat\t"
                           f"{pos['start'] + 1}\t{pos['end'] + 1}\t.\t"
                           f"{pos['orientation']}\t.\t"
                           f"ID={pos['kmer']};is_end={pos['is_end']}\n")
    
    def write_metadata(self, output_path):
        """Write k-mer statistics metadata"""
        with open(output_path, 'w') as f:
            f.write("kmer\ttotal_count\tend_count\n")
            end_counts = defaultdict(int)
            
            for positions in self.kmer_positions.values():
                for pos in positions:
                    if pos['is_end']:
                        end_counts[pos['kmer']] += 1
            
            for kmer, count in self.kmer_counts.items():
                f.write(f"{kmer}\t{count}\t{end_counts[kmer]}\n")

def main():
    parser = argparse.ArgumentParser(description='Find telomeric repeats in genome assemblies')
    parser.add_argument('input', help='Input FASTA file (can be gzipped)')
    parser.add_argument('--kmers', type=int, nargs='+', default=[6, 7],
                       help='K-mer sizes to analyze')
    parser.add_argument('--format', choices=['bed', 'gff'], default='bed',
                       help='Output format (default: bed)')
    parser.add_argument('--output-prefix', default='telomere_results',
                       help='Prefix for output files')
    
    args = parser.parse_args()
    
    # Initialize finder
    finder = TelomereFinder(args.kmers)
    
    # Read input file (handle both plain and gzipped FASTA)
    open_func = gzip.open if args.input.endswith('.gz') else open
    with open_func(args.input, 'rt') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            finder.process_sequence(record.id, record.seq)
    
    # Write outputs
    if args.format == 'bed':
        finder.write_bed(f"{args.output_prefix}.bed")
    else:
        finder.write_gff(f"{args.output_prefix}.gff")
    
    finder.write_metadata(f"{args.output_prefix}.stats.txt")
    
    # Display summary of top k-mers at ends
    print("\nTop telomeric repeats found at contig ends:")
    end_counts = defaultdict(int)
    for positions in finder.kmer_positions.values():
        for pos in positions:
            if pos['is_end']:
                end_counts[pos['kmer']] += 1
    
    for kmer, count in sorted(end_counts.items(), key=lambda x: x[1], reverse=True)[:5]:
        print(f"{kmer}: {count} occurrences at contig ends")

if __name__ == '__main__':
    main()
