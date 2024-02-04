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
