# this algorithm is not computationally efficient
# def fibonacci(n, k = 1):
#     if n == 0:
#         return 0
#     if n == 1:
#         return 1
#     else:
#         return (fibonacci((n - 1), k) + k * fibonacci((n - 2), k))


def fibonacci(n, k=1):
    '''
    Function to return the nth term of a fibonacci series
    Args:
        n: the value of n
        k: factor by which F(n-2) is multiplied with (no of rabbit pairs, one rabbits pair gives birth to)
    return:
        value: value of F(n)
    '''
    n_1, n_2 = 1, 1  # n_0 = 0
    for i in range(0, n):
        n_1, n_2 = n_2, (n_2 + (n_1 * k))
    return n_1


# def decaying_fibonacci(n, decay_rate, k=1):
#     '''
#     Function to return the nth term of a fibonacci series
#     Args:
#         n: the value of n
#         k: factor by which F(n-2) is multiplied with (no of rabbit pairs, one rabbits pair gives birth to)
#         decay_rate: number of terms a rabbit pair would be alive
#     return:
#         value: value of F(n)
#     '''

def fasta_to_dict(file_path):
    '''
    Function to convert text in a FASTA file to a pythin dictionery
    Args:
        file_path: dtring containing path to the fasta file
    Return:
        seq_dict: Dictionary for the fasta file sent, with keys as the sequence id and values as the sequence
    '''
    f1 = open(file_path, 'r')
    seq_dict = {}
    key = ""
    for line in f1:
        if ">" in line:
            # fasta header starts with ">" and ends with "\n"
            key = line[1:(len(line)-1)]
            seq_dict[key] = ""
        else:
            seq_dict[key] += line.strip('\n')
    return seq_dict


def motif_search_overlapping(seq, subseq, init_pos=0):
    '''
    Function to return the positions of occurances (1-indexed) of motifs in the genetic string
    Args:
        seq: the genetic string
        subseq: the motif
        init_pos: starting point in the seq
    Return:
        A list containing all the start_indices the motif occurs in the genetic string
    '''
    start_indices = []
    for i in range(init_pos, (len(seq) - len(subseq) + 1)):
        if seq[i:i + len(subseq)] == subseq:
            start_indices.append(i + 1)
    return start_indices


# pythonic way of motif_search_overlapping()
def CountOccurrences(string, substring):
    '''
    Function to return the number of occurances of substring in string"
    '''
    # Initialize count and start to 0
    count = 0
    start = 0

    # Search through the string till
    # we reach the end of it
    while start < len(string):

        # Check if a substring is present from
        # 'start' position till the end
        pos = string.find(substring, start)

        if pos != -1:
            # If a substring is present, move 'start' to
            # the next position from start of the substring
            start = pos + 1

            # Increment the count
            count += 1
        else:
            # If no further substring is present
            break
    # return the value of count
    return count


def most_frequent_kmer(seq, k):
    """
    Function to return the k-mer that is most frequently occuring in the sequence
    Args:
        seq: the genetic sequence
        k : length of k-mer
    Return:
        k-mers that occur the most times in the sequence
    """
    kmer_count = {}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if kmer in kmer_count.keys():
            kmer_count[kmer] += 1
        else:
            kmer_count[kmer] = 0
    max_value = max(kmer_count.items(), key=lambda x: x[1])[1]
    list_of_kmers = []
    for key, value in kmer_count.items():
        if value == max_value:
            list_of_kmers.append(key)

    return list_of_kmers


def clump_finder(seq, k, l, t):
    '''
    Function to identify k-mers that form clumps in the seq
    k-mer forms an (L, t)-clump in seq if there is an interval length L in which k-mer appears at least t times
    Args:
        seq: DNA sequence
        k: length of k-mer
        l: window size of clump
        t: frequency of kmer in window
    Return:
        List of all kmers that form (L,t) clump in seq
    '''
    kmer_clumps = []
    for i in range(len(seq) - l + 1):
        window = seq[i:i+l]
        for j in range(l - k + 1):
            current_kmer = window[j:j + k]
            count = CountOccurrences(window, current_kmer)
            if count >= t:
                kmer_clumps.append(current_kmer)
    return (set(kmer_clumps))


def min_skew_finder(seq):
    '''
    Function to return the positions (0 - indexed) of minimum skews in DNA seq
    Skew is the difference between the total number of occurrences of 'G' and 'C' in Genome.
    Args:
        Seq: DNA sequence
    Return:
        A list containing positions of occurances of minimum skew
    '''
    skew = []
    for i in range(len(seq) + 1):
        skew.append(seq[0:i].count("G") - seq[0:i].count("C"))
    return [i for i, x in enumerate(skew) if x == min(skew)]


def hamming_distance(seq_1, seq_2):
    '''Function to return the hamming distance between two DNA/RNA sequence'''
    ham_dist = 0
    if len(seq_1) != len(seq_2):
        print("Invalid Sequence sent. Make sure both sequences are of same length")
        return - 1
    else:
        for i in range(len(seq_1)):
            if seq_1[i] != seq_2[i]:
                ham_dist += 1
        return ham_dist


def motif_search_overlapping_approx(seq, subseq, d, init_pos=0):
    '''
    Function to return the positions of approximate occurances (0-indexed) of motifs in the genetic string
    Approximate occurance is when the motif occurs with less than d mutations. HammingDistance(occurance, motif) <= d
    Args:
        seq: the genetic string
        subseq: the motif
        d: tolerance for approximation
        init_pos: starting point in the seq
    Return:
        A list containing all the start_indices the motif occurs in the genetic string
    '''
    start_indices = []
    for i in range(init_pos, (len(seq) - len(subseq) + 1)):
        if (hamming_distance(seq[i:i + len(subseq)], subseq) <= d):
            start_indices.append(i)
    return start_indices
