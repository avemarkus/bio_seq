def filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0) -> dict:
    """
    Filters fastq sequences based on GC content, length, average quality.

    Args:
        seqs - dictionary containing fastq sequences.
        gc_bounds - GC content range. Default (0, 100), when all of the reads are included.
        length_bounds - length range. Default (0, 2**32).
        quality_threshold - average quality threshold, default 0 according to the phred33 scale.

    Returns:
        A dictionary with filtered sequences.
    """

    if not all(isinstance(value, tuple) and len(value) == 2 for value in seqs.values()):
        raise ValueError("Ошибка! Проверьте содержимое файла fastq.")

    filtered_seqs = {}
    
    for name, (sequence, quality) in seqs.items():
        gc_percentage = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
        if isinstance(gc_bounds, tuple):
            if gc_bounds[0] <= gc_percentage <= gc_bounds[1]:
                if isinstance(length_bounds, tuple):
                    if length_bounds[0] <= len(sequence) <= length_bounds[1]:
                        avg_quality = sum(ord(q) - 33 for q in quality) / len(quality)
                        if avg_quality >= quality_threshold:
                            filtered_seqs[name] = (sequence, quality)
                else:
                    if len(sequence) <= length_bounds:
                        avg_quality = sum(ord(q) - 33 for q in quality) / len(quality)
                        if avg_quality >= quality_threshold:
                            filtered_seqs[name] = (sequence, quality)
        else:
            if gc_percentage <= gc_bounds:
                if isinstance(length_bounds, tuple):
                    if length_bounds[0] <= len(sequence) <= length_bounds[1]:
                        avg_quality = sum(ord(q) - 33 for q in quality) / len(quality)
                        if avg_quality >= quality_threshold:
                            filtered_seqs[name] = (sequence, quality)
                else:
                    if len(sequence) <= length_bounds:
                        avg_quality = sum(ord(q) - 33 for q in quality) / len(quality)
                        if avg_quality >= quality_threshold:
                            filtered_seqs[name] = (sequence, quality)

    if not filtered_seqs:
        return("Удовлетворяющих последовательностей не найдено. Попробуйте другие аргументы.")
    
    return filtered_seqs