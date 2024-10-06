```markdown
# bio_seq Package

## Overview

The `bioseq` package is designed to facilitate bioinformatics operations on DNA and RNA sequences. It contains tools for filtering fastq files and performing sequence operations.

## Features

- **Filtering fastq sequences**
- **Transcription, reversal and complementation of sequences**
- **GC-persentage calculation**

## Usage

### Filtering FASTQ Sequences

```python
from bioseq.scripts.filter_fastq import filter_fastq

fastq_sequences = {
    'seq1': ('AGCTAGCTAGAT', 'IIIIIIIIIIII'),
    'seq2': ('GCGCGC', 'IIIIIIII')
}

filtered = filter_fastq(fastq_sequences, gc_bounds=(40, 60), length_bounds=(0, 100), quality_threshold=20)
print(filtered)
```

### DNA/RNA operations

```python
from bioseq.scripts.run_dna_rna_tools import run_dna_rna_tools

result = run_dna_rna_tools('aUGc', 'transcribe')
print(result)  # uACg
```

### Additional functions

You can also perform other operations such as `reverse`, `complement`, `reverse_complement`, and `count_GC` using the `run_dna_rna_tools` function.

```python
result = run_dna_rna_tools('ATGC', 'count_GC')
print(result)  # 50.0
```

## Conclusion

This package can be enhanced with additional functionalities and more advanced features as needed.

## P.S.

Guys, there were a lot of emergency situations. And they shot at me, and I fell into a puddle and remained there, because one good person gave us the task of creating this repository...

![:3](https://cs12.pikabu.ru/post_img/2020/10/29/11/160399913919458667.jpg)
