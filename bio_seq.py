"""

The bio_seq package is designed to facilitate bioinformatics operations on DNA and RNA sequences. 
It contains tools for filtering fastq files and performing sequence operations.

run_dna_rna_tools - takes nucleotide sequences and performs various operations on them.

filter_fastq - filters fastq sequences based on GC content, length, average quality.

"""

from .scripts_stock.run_dna_rna_tools import run_dna_rna_tools

from .scripts_stock.filter_fastq import filter_fastq
