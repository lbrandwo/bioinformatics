from itertools import product

"""
Bioinformatics Toolkit for Genomic Analysis

Introduction:

This Python script is a collection of bioinformatics functions tailored for genomic sequence analysis. It provides essential tools for pattern matching, 
mutation analysis, genomic feature identification, and sequence analysis. These functions facilitate a wide range of bioinformatics tasks, from basic 
sequence manipulation to advanced genomic inquiries, making it a valuable asset for researchers, educators, and developers in the field.

Capabilities:

1. **Pattern Matching**: Search for specific DNA sequences within a genome to locate genes, regulatory elements, and motifs.
2. **Mutation Analysis**: Compare DNA sequences to analyze genetic variations and mutations, aiding in evolutionary studies and disease research.
3. **Genomic Feature Identification**: Identify critical genomic features, such as replication origins, by analyzing nucleotide compositions and patterns.
4. **Sequence Analysis**: Perform detailed analyses of DNA sequences, including finding frequent k-mers with mismatches and computing reverse complements.

For examples and workflows check the /analysis folder.
"""


def PatternCount(Text, Pattern):
    """
    Counts the occurrences of a specific pattern within a given text.

    Parameters:
    Text (str): The text in which to search for the pattern.
    Pattern (str): The pattern to search for within the text.

    Returns:
    int: The number of times the pattern occurs in the text.
    """
    count = 0
    for i in range(0, (len(Text) - len(Pattern) + 1)):
        if Text[i: i + len(Pattern)] == Pattern:
            count += 1
    return count


def FrequencyTable(Text, k):
    """
    Generates a frequency map of all k-mers in a given text.

    Parameters:
    Text (str): The text from which to generate the frequency table.
    k (int): The length of the k-mers to search for.

    Returns:
    dict: A dictionary with k-mers as keys and their frequencies as values.
    """
    freqMap = {}
    n = len(Text)
    for i in range(n - k + 1):
        Pattern = Text[i:i+k]
        if Pattern not in freqMap:
            freqMap[Pattern] = 1
        else:
            freqMap[Pattern] += 1
    return freqMap


def MaxMap(freqMap):
    """
    Finds the maximum value in a frequency map.

    Parameters:
    freqMap (dict): A frequency map of k-mers and their counts.

    Returns:
    int: The maximum frequency found in the map.
    """
    return max(freqMap.values())


def BetterFrequentWords(Text, k):
    """
    Identifies the most frequent k-mers within a given text, including those with the same maximum frequency.

    Parameters:
    Text (str): The text to analyze.
    k (int): The length of the k-mers.

    Returns:
    str: A space-separated string of the most frequent k-mers.
    """
    FrequentPatterns = []
    freqMap = FrequencyTable(Text, k)
    maxVal = MaxMap(freqMap)
    for Pattern in freqMap:
        if freqMap[Pattern] == maxVal:
            FrequentPatterns.append(Pattern)
    return " ".join(FrequentPatterns)


def reverse_complement(pattern):
    """
    Computes the reverse complement of a DNA sequence.

    Parameters:
    pattern (str): A DNA sequence.

    Returns:
    str: The reverse complement of the given DNA sequence.
    """
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_comp = ''.join([complement_map[nuc] for nuc in reversed(pattern)])
    return reverse_comp


def PatternMatcher(Pattern, Genome):
    """
    Finds all starting positions of a pattern within a given genome sequence.

    Parameters:
    Pattern (str): The DNA pattern to search for.
    Genome (str): The genome sequence to search within.

    Returns:
    str: A space-separated string of starting positions where the pattern is found.
    """
    output = []
    for i in range(len(Genome) - len(Pattern) + 1):
        if Genome[i: i + len(Pattern)] == Pattern:
            output.append(str(i))
    return " ".join(output)


def find_substring_positions(genome, pattern):
    """
    Identifies the start positions where a pattern is found within a given genome.

    Parameters:
    genome (str): The genome to search.
    pattern (str): The pattern to find within the genome.

    Returns:
    str: A space-separated string of start positions where the pattern occurs in the genome.
    """
    positions = []
    for i in range(len(genome) - len(pattern) + 1):
        if genome[i:i+len(pattern)] == pattern:
            positions.append(str(i))
    return ' '.join(positions)


def FindClumps(Text, k, L, t):
    """
    Finds patterns forming clumps within a window of a genome.

    Parameters:
    Text (str): The genome sequence.
    k (int): Length of k-mer.
    L (int): Length of the window to examine for clumps.
    t (int): Minimum number of occurrences to qualify as a clump.

    Returns:
    str: A space-separated string of unique k-mers that form clumps within the window.
    """
    Patterns = set()
    n = len(Text)
    for i in range(n - L + 1):
        Window = Text[i:i+L]
        freqMap = FrequencyTable(Window, k)
        for s in freqMap:
            if freqMap[s] >= t:
                Patterns.add(s)
    return ' '.join(Patterns)


def compute_skew(genome):
    """
    Computes the skew (#G - #C) for every position in the genome.

    Parameters:
    genome (str): The genome sequence.

    Returns:
    list: Skew values for each position in the genome.
    """
    skew = [0]
    for i in range(len(genome)):
        if genome[i] == 'G':
            skew.append(skew[-1] + 1)
        elif genome[i] == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew


def find_minimum_skew_positions(genome):
    """
    Finds positions in a genome where the skew diagram reaches a minimum.

    Parameters:
    genome (str): The genome sequence.

    Returns:
    list: Positions where the skew reaches its minimum value.
    """
    skew_values = compute_skew(genome)
    min_skew = min(skew_values)
    return [i for i, value in enumerate(skew_values) if value == min_skew]


def HammingDistance(p, q):
    """
    Calculates the Hamming distance between two strings.

    Parameters:
    p (str): First string.
    q (str): Second string.

    Returns:
    int: The Hamming distance between the two strings.
    """
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1
    return count


def find_substring_positions(genome, pattern, d):
    """
    Finds all positions where a pattern occurs in a genome, allowing for up to d mismatches.

    Parameters:
    genome (str): The genome to search.
    pattern (str): The pattern to find within the genome.
    d (int): Maximum number of allowed mismatches.

    Returns:
    int: The number of positions where the pattern occurs with up to d mismatches.
    """
    positions = []
    pattern_length = len(pattern)
    for i in range(len(genome) - pattern_length + 1):
        if HammingDistance(genome[i:i+pattern_length], pattern) <= d:
            positions.append(str(i))
    return len(positions)


def hamming_distance(s1, s2):
    """
    Calculates the Hamming distance between two DNA sequences.

    Parameters:
    s1 (str): The first DNA sequence.
    s2 (str): The second DNA sequence.

    Returns:
    int: The number of mismatches between the two sequences.
    """
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def mismatches(kmer, d):
    """
    Generates all possible k-mers that are within a Hamming distance of d from the original k-mer.

    Parameters:
    kmer (str): The original k-mer.
    d (int): The Hamming distance.

    Returns:
    generator: A generator yielding all k-mers within the specified distance.
    """
    nucleotides = 'ACGT'
    for positions in product(nucleotides, repeat=len(kmer)):
        if hamming_distance(kmer, ''.join(positions)) <= d:
            yield ''.join(positions)


def frequent_words_with_mismatches(text, k, d):
    """
    Identifies the most frequent k-mers within a text, allowing for up to d mismatches.

    Parameters:
    text (str): The DNA sequence.
    k (int): The length of the k-mers.
    d (int): The allowed number of mismatches.

    Returns:
    list: A list of the most frequent k-mers within the specified mismatch tolerance.
    """
    possible_kmers = set()
    kmer_counts = {}

    for i in range(len(text) - k + 1):
        window = text[i:i+k]
        for variant in mismatches(window, d):
            possible_kmers.add(variant)

    for kmer in possible_kmers:
        kmer_counts[kmer] = sum(hamming_distance(
            text[i:i+k], kmer) <= d for i in range(len(text) - k + 1))

    max_count = max(kmer_counts.values())
    return [kmer for kmer, count in kmer_counts.items() if count == max_count]


def Neighbors(pattern, d):
    """
    Generates all k-mers that differ from the given k-mer by at most d mismatches.

    Parameters:
    pattern (str): The original k-mer.
    d (int): The maximum number of allowed mismatches.

    Returns:
    set: A set of all possible k-mers within the specified Hamming distance.
    """
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    neighbors = set()
    suffix_neighbors = Neighbors(pattern[1:], d)
    for text in suffix_neighbors:
        if HammingDistance(pattern[1:], text) < d:
            for nucleotide in 'ACGT':
                neighbors.add(nucleotide + text)
        else:
            neighbors.add(pattern[0] + text)
    return neighbors


def FrequentWordsWithMismatchesAndReverse(genome, k, d):
    """
    Finds the most frequent k-mers with up to d mismatches in a genome, including the reverse complement.

    Parameters:
    genome (str): The genome sequence.
    k (int): The k-mer length.
    d (int): The allowed number of mismatches.

    Returns:
    list: The most frequent k-mers with mismatches and their reverse complements.
    """
    freqMap = {}
    for i in range(len(genome) - k + 1):
        pattern = genome[i:i+k]
        neighborhood = Neighbors(pattern, d) | Neighbors(
            ReverseComplement(pattern), d)
        for neighbor in neighborhood:
            freqMap[neighbor] = freqMap.get(neighbor, 0) + 1

    maxCount = max(freqMap.values())
    return [kmer for kmer, count in freqMap.items() if count == maxCount]


def compute_skew(genome):
    """
    Computes the skew array for a genome, indicating the difference between the cumulative number of 'G' and 'C'.

    Parameters:
    genome (str): The genome sequence.

    Returns:
    list: The skew values for each position in the genome.
    """
    skew = [0]  # Skew0 is 0
    for i in range(len(genome)):
        if genome[i] == 'G':
            skew.append(skew[-1] + 1)
        elif genome[i] == 'C':
            skew.append(skew[-1] - 1)
        else:
            skew.append(skew[-1])
    return skew


def find_minimum_skew_positions(genome):
    """
    Identifies positions in the genome where the skew reaches a minimum, potential oriC locations.

    Parameters:
    genome (str): The genome sequence.

    Returns:
    list: Positions of minimum skew in the genome.
    """
    skew_values = compute_skew(genome)
    min_skew = min(skew_values)
    return [i for i, value in enumerate(skew_values) if value == min_skew]
