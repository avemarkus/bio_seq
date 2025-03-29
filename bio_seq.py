from abc import ABC, abstractmethod
from typing import Dict

class BiologicalSequence(ABC):
    def __init__(self, sequence: str):
        self.sequence = sequence
        self.validate_alphabet()  # сразу проверяем алфавит

    @abstractmethod
    def validate_alphabet(self):
        """check if sequence is ok"""
        pass

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.sequence}')"

class NucleicAcidSequence(BiologicalSequence):
    @abstractmethod
    def _complement_base(self, base):
        """get complement base"""
        pass

    def complement(self):
        complemented = ''.join(self._complement_base(b) for b in self.sequence)
        return self.__class__(complemented) # возвращаем тот же класс

    def reverse(self):
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self):
        return self.complement().reverse()

    def count_gc(self):
        g_count = self.sequence.upper().count('G')
        c_count = self.sequence.upper().count('C')
        length = len(self.sequence)
        if length > 0:
            return (g_count + c_count) / length * 100
        return 0.0

class DNASequence(NucleicAcidSequence):
    def validate_alphabet(self):
        valid = {'A', 'T', 'G', 'C'}  # только эти буквы ок
        for base in self.sequence:
            if base not in valid:
                raise ValueError(f"Oops, '{base}' isn't a valid DNA base!")

    def _complement_base(self, base):
        complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return complements[base]

    def transcribe(self):
        rna_seq = self.sequence.replace('T', 'U') # заменяем T на U
        return RNASequence(rna_seq)

class RNASequence(NucleicAcidSequence):
    def validate_alphabet(self):
        valid = {'A', 'U', 'G', 'C'}
        for x in self.sequence:
            if x not in valid:
                raise ValueError(f"Hey, '{x}' doesn't belong in RNA!")

    def _complement_base(self, base):
        comp = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}  # комплементарные пары
        return comp[base]

class AminoAcidSequence(BiologicalSequence):
    def validate_alphabet(self):
        """Sequence contains only standard amino acids?"""
        amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
        for aa in self.sequence:
            if aa not in amino_acids:
                raise ValueError(f"'{aa}' is not a valid amino acid!")

    def calculate_hydrophobicity(self):  # считаем гидрофобность
        hydrophobicity_dict = { 'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
                               'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
                               'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
                               'S': -0.8 ,'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3}
        total = sum(hydrophobicity_dict[x] for x in self.sequence)
        return total / len(self.sequence)