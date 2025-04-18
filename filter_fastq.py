from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
from typing import Dict, Tuple

def filter_fastq(fastq_file: str, gc_bounds: Tuple[float, float] = (0, 100), 
                 length_bounds: Tuple[int, int] = (0, 2**32), quality_threshold: float = 0) -> Dict[str, Tuple[str, str]]:
    """
    Filter FastQ sequences based on GC content, length, and average quality.

    Args:
        fastq_file: Path to the FastQ file.
        gc_bounds: GC content range in percentage. Defaults to (0, 100).
        length_bounds: Sequence length range. Defaults to (0, 2**32).
        quality_threshold: Minimum average quality threshold (Phred scale). Defaults to 0.

    Returns:
        A dictionary with filtered sequences in the format {name: (sequence, quality)}.
        Returns a string message if no sequences meet the criteria.
    """
    filtered_seqs: Dict[str, Tuple[str, str]] = {}

    try:
        with open(fastq_file, 'r'):
            pass
    except FileNotFoundError:
        raise FileNotFoundError(f"FastQ file '{fastq_file}' not found")

    for record in SeqIO.parse(fastq_file, "fastq"):
        sequence = str(record.seq)
        gc_content = GC(record.seq)
        seq_length = len(sequence)
        quality_scores = record.letter_annotations["phred_quality"]
        avg_quality = np.mean(quality_scores) if quality_scores else 0.0

        if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
            length_bounds[0] <= seq_length <= length_bounds[1] and
            avg_quality >= quality_threshold):
            quality_str = "".join(chr(q + 33) for q in quality_scores)
            filtered_seqs[record.id] = (sequence, quality_str)

    if not filtered_seqs:
        return "Удовлетворяющих последовательностей не найдено. Попробуйте другие аргументы."
    
    return filtered_seqs
