# from bioinformatics_tools import *


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


def edge_list_format_reader(file_path):
    '''
    Function retrieve nodes and edges from file in edge list format
    Args:
        file_path: a string containing the file path
    Return:
        n: number of nodes
        e: number of edges
        edges: a list containing all edges(2 nodes connected by a edge)
    '''
    n = 0
    e = 0
    edges = []
    f = open(file_path, "r")
    start = -1
    for i in f:
        start = start + 1
        if start == 0:
            graph_params = i.strip().split()
            n = int(graph_params[0])
            e = int(graph_params[1])
        if start > 0:
            char_edge = i.strip().split()
            int_edge = [int(j) for j in char_edge]
            edges.append(int_edge)
    return n, e, edges


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


def kmer_frequencies(seq, k):
    '''
    Function to return the frequency of occurances of all possible k-mers
    Args:
        seq: the DNA sequence
        k: length of k-mer
    Return:
        A list of frequencies of all kmer ordered lexicographically
    '''
    from itertools import product
    all_kmers = list(product(["A", "C", "G", "T"], repeat=k))
    frequencies = []
    for i in all_kmers:
        kmer = "".join(i)
        frequencies.append(CountOccurrences(seq, kmer))
    return frequencies


def lexicographic_kmer(index, k):
    '''
    Function to return the kmer that has a rank = index when ordered lexicographically
    Args:
        index: rank in lexicographic order
        k: length of kmer
    Return:
        the kmer of rank = index
    '''
    ranked_nucleotide = {0: "A", 1: "C", 2: "G", 3: "T"}
    if k == 1:
        return ranked_nucleotide[index]
    return lexicographic_kmer((index//4), (k-1)) + ranked_nucleotide[index % 4]


def lexicographic_kmer_rank(kmer):
    '''
    Function to return the rank of kmer in the lexicographic order
    Arg:  
        kmer: the kmer
    Return:
        rank of given kmer among all kmers ordered lexicographically
    '''
    nucleotide_rank = {"A": 0, "C": 1, "G": 2, "T": 3}
    if len(kmer) == 1:
        return nucleotide_rank[kmer]
    return 4*lexicographic_kmer_rank(kmer[:-1]) + nucleotide_rank[kmer[-1]]


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
            kmer_count[kmer] = 1
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


def kmer_neighbours(kmer, d):
    '''
    Function to return the kmers that has a hamming distance less than d from the given kmer
    Args:
        kmer: the kmer around which the neighbourhood needs to be calculated
        d: the size of the neighbourhood
    Returns:
        A list of kmers in the d-neighbourhood
    '''
    if d == 0:
        return [kmer]
    if len(kmer) == 1:
        return ["A", "C", "G", "T"]
    neighbourhood = []
    suffix_neighbors = kmer_neighbours(kmer[1:], d)
    for i in suffix_neighbors:
        if hamming_distance(i, kmer[1:]) < d:
            for j in ["A", "C", "G", "T"]:
                neighbourhood.append((j + i))
        else:
            neighbourhood.append((kmer[0] + i))
    return neighbourhood


def most_frequent_kmer_approx(seq, k, d):
    """
    Function the return the kmers with atmost d-mismatches in the string
    Args:
        seq: the DNA sequence
        k: length of kmer
        d: maximum mismatches in kmer allowed
    Return:
        The list of kmers that approximately occur in the seq
    """
    kmer_count = {}
    kmer_neighbourhood = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmer_neighbourhood += kmer_neighbours(kmer, d)
    for kmer in kmer_neighbourhood:
        if kmer in kmer_count.keys():
            kmer_count[kmer] += 1
        else:
            kmer_count[kmer] = 1
    max_value = max(kmer_count.items(), key=lambda x: x[1])[1]
    list_of_kmers = []
    for key, value in kmer_count.items():
        if value == max_value:
            list_of_kmers.append(key)

    return list_of_kmers


def most_frequent_kmer_approx_reverse(seq, k, d):
    """
    Function the return the kmers and its reverse compliment with atmost d-mismatches in the string
    Args:
        seq: the DNA sequence
        k: length of kmer
        d: maximum mismatches in kmer allowed
    Return:
        The list of kmers that approximately occur in the seq
    """
    kmer_count = {}
    kmer_neighbourhood = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmer_neighbourhood += kmer_neighbours(kmer, d)
        kmer_reverse = dna_reverse_compliment(kmer)
        kmer_neighbourhood += kmer_neighbours(kmer_reverse, d)
    for kmer in kmer_neighbourhood:
        if kmer in kmer_count.keys():
            kmer_count[kmer] += 1
        else:
            kmer_count[kmer] = 1
    max_value = max(kmer_count.items(), key=lambda x: x[1])[1]
    list_of_kmers = []
    for key, value in kmer_count.items():
        if value == max_value:
            list_of_kmers.append(key)

    return list_of_kmers


def motif_random_probability(motif, gc_content):
    '''
    Function to return the probability of a motif occuring in random given the gc_content of the complete sequence
    If the GC-content is x, then we set the symbol probability of C and G equal to x/2 
    and the symbol probability of A and T equal to (1−x)/2
    Args:
        seq: the motif
        gc_content: GC content of the sequence
    Returns:
        The probability that the motif occured by random chance
    '''
    symbol_probability = {"A": (1 - gc_content) / 2,
                          "C": gc_content / 2,
                          "T": (1 - gc_content) / 2,
                          "G": gc_content / 2}
    probability = 1
    print(symbol_probability)
    for i in motif:
        probability = probability * symbol_probability[i]
    return probability


def simple_graph_adjacency_matrix(n, e, edges):
    '''
    Function to return adjecency matrix of a simple graph
    Args:
        n: Number of nodes
        e: Number of edges
        edges: the undirectional vertices
    Return:
        The graph in a n*n matrix format (entries of the matrix mention connectivity b/w nodes)
    '''
    adjacency_matrix = [[0 for i in range(n)] for j in range(n)]

    for i in range(e):
        m = edges[i][0] - 1
        n = edges[i][1] - 1
        adjacency_matrix[m][n] = 1
        adjacency_matrix[n][m] = 1

    return adjacency_matrix


def simple_graph_adjacency_list(n, e, edges):
    '''
    Function to return adjecency matrix of a simple graph
    Args:
        n: Number of nodes
        e: Number of edges
        edges: the undirectional vertices
    Return:
        The graph in a n*n matrix format (entries of the matrix mention connectivity b/w nodes)
    '''
    adjacency_list = {}

    # initializing adjacency list
    for i in range(n):
        adjacency_list[i + 1] = []

    for vertex in edges:
        node_i = vertex[0]
        node_j = vertex[1]
        adjacency_list[node_i].append(node_j)
        adjacency_list[node_j].append(node_i)

    return adjacency_list


def simple_graph_double_degree(adjacency_matrix, adjacency_list):
    '''
    Function to return double degree of all nodes in a simple graph
    Args:
        adjacency_matrix: the adjacency matrix of the simple graph
        adjacency_list: the adjacency list of the simple graph
    Return:
        list containing double degrees of all nodes
    '''
    double_degrees = [0 for i in range(len(adjacency_matrix))]

    for i in range(len(double_degrees)):
        neighbouring_nodes = adjacency_list[i + 1]
        for node in neighbouring_nodes:
            double_degrees[i] += sum(adjacency_matrix[node - 1])

    return double_degrees
