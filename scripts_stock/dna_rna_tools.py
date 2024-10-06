def run_dna_rna_tools(*arg: str) -> str or float or list:
    """
    Function takes nucleotide sequences and performs
    various operations on them.

    Args:
        *arg: Nucleotide sequences and the name of the operation (optionally) on the last position.
        The list of operations:
        transcribe — return the transcribed sequence
        reverse — return the reversed sequence
        complement - return the complementary sequence
        reverse_complement — return the reverse complementary sequence
        count_GC - calculates the GC-percentage in a sequence
        
    Returns:
        The result of performing an operation.
        If no operation is specified, returns input sequences.
    """
    if len(arg) < 1:
        return "Введите аргументы функции!"

    if len(arg) > 1 and arg[-1] in ['transcribe', 'reverse', 'complement', 'reverse_complement', 'count_GC']:
        procedure = arg[-1]
        sequences = arg[:-1]
    else:
        procedure = None
        sequences = arg

    for sequence in sequences:
        if not is_norm_sequence(sequence):
            raise ValueError("Некорректная последовательность: " + sequence)

    results = []

    for sequence in sequences:
        if procedure == 'transcribe':
            results.append(transcribe(sequence))
        elif procedure == 'reverse':
            results.append(reverse(sequence))
        elif procedure == 'complement':
            results.append(complement(sequence))
        elif procedure == 'reverse_complement':
            results.append(reverse_complement(sequence))
        elif procedure == 'count_GC':
            results.append(count_GC(sequence))
        else:
            results = list(sequences)

    if len(results) == 1:
        return results[0]
    return results


def transcribe(sequence: str) -> str:
    """
    Return the transcribed sequence.

    Args:
        Initial sequence as str.

    Returns:
        Transcribed sequence as str
    """
    return sequence.replace("T", "U").replace("t", "u")


def reverse(sequence: str) -> str:
    """
    Return the reversed sequence.

    Args:
        Initial sequence as str.

    Returns:
        Reversed sequence as str
    """
    return sequence[::-1]


def complement(sequence: str) -> str:
    """
    Return the complementary sequence.

    Args:
        Initial sequence as str.

    Returns:
        Complementary sequence as str.
    """
    complement_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "U": "A",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
        "u": "a",
    }
    return "".join(complement_dict[nucleotide] for nucleotide in sequence)


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complementary sequence.

    Args:
        Initial sequence as str.

    Returns:
        Complementary and reverse sequence as str.
    """
    return reverse(complement(sequence))


def count_GC(sequence: str) -> float:
    """
    Calculates the GC-percentage in a sequence.

    Args:
        Initial sequence as str.

    Returns:
        GC-percentage in a sequence as float.
    """
    g_count = sequence.upper().count("G")
    c_count = sequence.upper().count("C")
    total_count = len(sequence)

    if total_count == 0:
        return 0

    gc_percentage = (g_count + c_count) / total_count * 100
    return gc_percentage


def is_norm_sequence(sequence: str) -> bool:
    """
    Checks whether the sequence is correct.

    Args:
        Initial sequence as str.

    Returns:
        True if the sequence is correct DNA or RNA, False otherwise.
    """
    dna = set("ATCGatcg")
    rna = set("AUCGucgu")

    return set(sequence).issubset(dna) or set(sequence).issubset(rna)
