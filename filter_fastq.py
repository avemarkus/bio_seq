import argparse
import logging
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import numpy as np
from typing import Dict, Tuple

logger = logging.getLogger(__name__) # логгирование
logger.setLevel(logging.INFO)
if not logger.handlers:
    file_handler = logging.FileHandler("filter.log")
    file_handler.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    
def filter_fastq(fastq_file: str, gc_bounds: Tuple[float, float] = (0, 100),
                 length_bounds: Tuple[int, int] = (0, 2**32), quality_threshold: float = 0) -> Dict[str, Tuple[str, str]]:
    """
    Filter sequences based on GC content, length, and average quality.

    Args:
        fastq_file: Path to the fastq file.
        gc_bounds: GC content range in percentage (min-max). Defaults to (0, 100).
        length_bounds: Sequence length range (min-max). Defaults to (0, 2**32).
        quality_threshold: Minimum average quality threshold (Phred scale). Defaults to 0.

    Returns:
        A dictionary with filtered sequences in the format {name: (sequence, quality)}.
        Returns a string message if no sequences meet the criteria or if the file is invalid.

    Raises:
        FileNotFoundError: If the fastq file does not exist.
    """
    filtered_seqs: Dict[str, Tuple[str, str]] = {}

    logger.info(f"Starting filtering of fastq file: {fastq_file}")
    logger.info(f"Parameters: GC={gc_bounds}, Length={length_bounds}, Quality={quality_threshold}")

    try:
        with open(fastq_file, 'r'):
            pass
    except FileNotFoundError:
        logger.error(f"Fastq file '{fastq_file}' was not found")
        raise FileNotFoundError(f"Fastq file '{fastq_file}' was not found")

    try:
        for record in SeqIO.parse(fastq_file, "fastq"):
            sequence = str(record.seq)
            gc_content = gc_fraction(record.seq) * 100
            seq_length = len(sequence)
            quality_scores = record.letter_annotations["phred_quality"]
            avg_quality = np.mean(quality_scores) if quality_scores else 0.0

            if (gc_bounds[0] <= gc_content <= gc_bounds[1] and
                length_bounds[0] <= seq_length <= length_bounds[1] and
                avg_quality >= quality_threshold):
                quality_str = "".join(chr(q + 33) for q in quality_scores)
                filtered_seqs[record.id] = (sequence, quality_str)
    except ValueError as e:
        logger.error(f"Invalid format in {fastq_file}: {str(e)}")
        return "Invalid format. Please check the file."

    if not filtered_seqs:
        logger.info("No sequences was found the filtering criteria.")
        return "No sequences was found the filtering criteria. Try another parameters."

    logger.info(f"Filtered {len(filtered_seqs)} sequences")
    return filtered_seqs

def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Filter fastq sequences by GC content, length, and quality.")
    parser.add_argument("fastq_file", type=str, help="Path to the input fastq file")
    parser.add_argument("--gc-bounds", type=float, nargs=2, default=(0, 100),
                        help="GC content range (min-max, percentage), default: 0-100")
    parser.add_argument("--length-bounds", type=int, nargs=2, default=(0, 2**32),
                        help="Sequence length range (min-max), default: 0-2^32")
    parser.add_argument("--quality-threshold", type=float, default=0,
                        help="Minimum average Phred quality score, default: 0")
    return parser.parse_args()

def main():
    """
    Main function to run the fastq filter.
    """
    args = parse_arguments()
    result = filter_fastq(
        fastq_file=args.fastq_file,
        gc_bounds=tuple(args.gc_bounds),
        length_bounds=tuple(args.length_bounds),
        quality_threshold=args.quality_threshold
    )
    if isinstance(result, str):
        print(result)
    else:
        print(f"{len(result)} sequences were filtered:")
        for seq_id, (seq, qual) in result.items():
            print(f">{seq_id}\n{seq}\nQuality: {qual}")

if __name__ == "__main__":
    main()