# telomere repo
Scripts to identify putative telomeric and subtelomeric sequences from draft genome assemblies.

# Telomere finder tutorial

This tool helps identify potential telomeric repeats in genome assemblies by analyzing k-mer occurrences at contig ends. It supports both compressed and uncompressed FASTA inputs and can output results in BED or GFF format.

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/telomere-finder
cd telomere-finder

# Install dependencies
pip install biopython
```

## Usage

Basic usage with default settings (6-mers and 7-mers, BED output):

```bash
python telomere_finder.py input.fasta
```

Analyze specific k-mer sizes and output GFF:

```bash
python telomere_finder.py --kmers 6 7 8 --format gff input.fasta
```

## Demo Example

Let's analyze a simple demo genome with known telomeric repeats. The demo file `demo.fasta` contains three chromosomes:
- chr1: Contains multiple TTTAGGG repeats at both ends
- chr2: Has TTAGGG repeats internally and CCCTAA at the end
- chr3_reverse: Contains the canonical telomere sequence in reverse orientation (CCCTAA)

```bash
# Run analysis on demo genome
python telomere_finder.py --kmers 6 7 demo.fasta
```

Expected output:
```
Top telomeric repeats found at contig ends:
TTAGGG: 12 occurrences at contig ends
TTTAGGG: 8 occurrences at contig ends
CCCTAA: 12 occurrences at contig ends
```

## Output Files

1. `.bed` or `.gff` file containing positions of all k-mer occurrences:
   - BED format (0-based): `chromosome start end kmer score strand`
   - GFF format (1-based): Includes additional metadata about end proximity

2. `.stats.txt` file containing:
   - K-mer occurrence counts
   - Number of times each k-mer appears near contig ends
   - Distribution of k-mers across the genome

## Finding Telomeres

The tool is particularly useful for:
1. Identifying potential telomeric repeats in draft assemblies
2. Validating genome assembly completeness
3. Detecting non-canonical telomere sequences
4. Analyzing telomere orientation and distribution

Note that the tool considers sequences within 1000bp of contig ends as potential telomeric regions by default.

## Advanced Usage

### Compressed Input Files

The tool automatically handles gzipped FASTA files:

```bash
python telomere_finder.py input.fasta.gz
```

### Custom End Region Size

Modify the `end_region_size` parameter in the code to adjust what's considered an "end region" (default: 1000bp).

### Output Format Conventions

- BED format: 0-based, half-open intervals
- GFF format: 1-based, closed intervals

## Common Telomeric Sequences

Human and many vertebrates:
- TTAGGG (canonical)
- CCCTAA (reverse complement)

Some known variants:
- TTTAGGG
- TTAGGGG
- TTAGG

## Tips for Analysis

1. Start with multiple k-mer sizes to catch variants
2. Look for reverse complement sequences
3. Pay attention to k-mer clustering at contig ends
4. Consider both frequency and position of repeats

## Contributing

Contributions are welcome! Please submit issues and pull requests to the GitHub repository.
