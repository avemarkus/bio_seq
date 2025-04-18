# bio_seq Package

## Intro

The `bioseq` package is designed to facilitate bioinformatics operations on DNA, RNA and amino acid sequences and provides a robust object-oriented framework for sequence manipulation.
 
## Features

- **Filtering fastq sequences**
- **Operations on DNA and RNA sequences: transcription, complementation, reversal, reverse complementation, and GC content calculation**
- **Analysis of aminoacid sequences, including hydrophobicity calculationn**


Before usage ensure you have Python installed with the required dependencies (`biopython`, `numpy`):
```bash
pip install biopython numpy
```

## Usage

### Filtering FASTQ Sequences

```python

filtered = filter_fastq("example.fastq", gc_bounds=(40, 60), length_bounds=(50, 150), quality_threshold=20)
if isinstance(filtered, dict):
    for name, (seq, qual) in filtered.items():
        print(f"ID: {name}, Seq: {seq[:10]}..., Quality: {qual[:10]}...")
else:
    print(filtered)
    
```

### DNA/RNA operations

```python
from bio_seq import DNASequence, RNASequence

# DNA example
dna = DNASequence("GATTACA")
print(dna)                    # GATTACA
print(dna.complement())       # CTAATGT
print(dna.reverse())          # ACATTAG
print(dna.reverse_complement()) # TGTAATC
print(dna.count_gc())         # 28.57
rna = dna.transcribe()
print(rna)                    # GAUUACA

# RNA example
rna = RNASequence("AUGCGU")
print(rna)                    # AUGCGU
print(rna.complement())       # UACGCU
```

### Aminoacid operations

```python
from bio_seq import AminoAcidSequence

protein = AminoAcidSequence("MILVFW")
print(protein)                    # MILVFW
print(protein.calculate_hydrophobicity())  # 2.72
```

## Conclusion

This bio_seq package now based an OOP approach, making bioinformatics analyses more modular and extensible nevertheless it can be enhanced with additional functionalities and more advanced features as needed.

## P.S.S.

Guys, I can't take it anymore...

![:3](https://www.meme-arsenal.com/memes/ff70fc4d019a6fda321cf501f7b2ab32.jpg)
