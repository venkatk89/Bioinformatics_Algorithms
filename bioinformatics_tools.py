# importing structures in bio_struct file
from bio_structures import *
from algorithmic_tools import *
import numpy as np


# Useful functions in bioinformatics


def dna_seq_validate(seq):
    '''
    Function to validate a seq
    Args: The DNA sequence
    Return: The processed string if valid DNA sequence
            0 if invalid DNA sequence

    '''
    temp_seq = seq.upper()
    for nucleotide in temp_seq:
        if nucleotide not in Nucleotides:
            return False
    return temp_seq


def dna_freq_counter(seq):
    '''
    Function to find the length of a DNA or RNA sequence
    Args: The DNA sequence
    Return: dict containing the count of each nucleotides
    '''
    temp_seq = seq.upper()
    temp_counter = {"A": 0, "T": 0, "C": 0, "G": 0}
    for nucleotide in temp_seq:
        temp_counter[nucleotide] = temp_counter[nucleotide] + 1
    return temp_counter


def rna_transcription(seq):
    '''
    Function to transcribe DNA to RNA
    Args: The DNA sequence
    Return: The RNA sequence transcribed from the DNA
    '''
    temp_seq = seq.upper()
    return temp_seq.replace("T", "U")


def dna_reverse_compliment(seq):
    '''
    Function to produce the reverse complement of a dna strand
    Args: The DNA sequence
    Return: The reverse compliment of the DNA strand. (BP is complimented, then string is reversed)
    '''
    temp_seq = ""
    for i in range(len(seq)):
        temp_seq += Base_pairs[seq[(len(seq) - 1 - i)]]
    return temp_seq


def gc_content(seq):
    '''
    Function to return the GC count of a DNA/RNA sequence
    Args:
        seq: A valid DNA/RNA sequence
    Return:
        gc: The gc count of seq
    '''
    return round(((seq.count("C") + seq.count("G")) / len(seq)), 6)


def gc_content_subset(seq, k):
    '''
    Function to return the GC content of substrings in a DNA/RNA sequence
    Ags:
        Seq: A valid DNA/RNA sequence
        k: length of substring. It'll be the window size
    Return:
        res: a list containing gc content of each windows
    '''
    res = []
    for i in range(0, (len(seq) - k + 1), k):
        substring = seq[i:i+k]
        res.append(gc_content(substring))
    return res


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


def dna_translation(seq, init_pos=0):
    '''
    Function to Translate DNA strand to corresponding Protein string
    Args:
        seq: DNA sequence
        init_pos: position from which translation should be started
    Return:
        A Protein string
    '''
    protein_string = ""
    for i in range(init_pos, len(seq) - 2, 3):
        protein_string += DNA_Codons[seq[i:i + 3]]
    return protein_string


def rna_translation(seq, init_pos=0):
    '''
    Function to Translate RNA strand to corresponding Protein string
    Args:
        seq: RNA sequence
        init_pos: position from which translation should be started
    Return:
        A Protein string
    '''
    protein_string = ""
    for i in range(init_pos, len(seq) - 2, 3):
        protein_string += RNA_Codons[seq[i:i + 3]]
    return protein_string


def generate_reading_frames(seq):
    '''
    Function to translate all 6 combos of the sequence (3 in seq + 3 in reverse complement of seq)
    Args:
        seq: A valid DNA string
    Return:
        A list of 6 protein strings
    '''
    frames = []
    frames.append(dna_translation(seq, 0))
    frames.append(dna_translation(seq, 1))
    frames.append(dna_translation(seq, 2))
    frames.append(dna_translation(dna_reverse_compliment(seq), 0))
    frames.append(dna_translation(dna_reverse_compliment(seq), 1))
    frames.append(dna_translation(dna_reverse_compliment(seq), 2))
    return frames


def proteins_from_seq(aa_seq):
    '''
    Function to return possible proteins from an amino acid sequence
    Args:
        aa_seq: An amino acid sequence
    Return:
        A list of all possible proteins
    '''
    all_proteins = []
    current_protein_tracking = []
    for aa in list(aa_seq):
        if aa == "_":   # if stop codon is found
            for p in current_protein_tracking:
                all_proteins.append(p)
            current_protein_tracking = []

        else:
            if aa == "M":
                current_protein_tracking.append("")
            for i in range(len(current_protein_tracking)):
                current_protein_tracking[i] += aa
    return all_proteins


def proteins_from_seq_orfs(seq, init_pos=0, end_pos=3, ordered=False):
    '''
    Function to find all possible proteins from open reading frames of a sequence
    Args:
        seq: A valid DNA sequence
        init_pos: starting position to consider in the seq
        end_pos: ending position to consider in the seq
        ordered: whether to sort the list of all proteins
    Return:
        A list of all possible proteins from open reading frames of a sequence
    '''
    all_proteins = []
    if init_pos >= end_pos:
        print("Enter valid start and end positions")
    else:
        # generate all 6 open reading frames from sequence
        reading_frames = generate_reading_frames(seq[init_pos:end_pos])
        # find possible proteins from each reading frame and add it to the list
        for frame in reading_frames:
            proteins = proteins_from_seq(frame)
            for each_protein in proteins:
                all_proteins.append(each_protein)
        if ordered == True:
            return sorted(all_proteins, key=len, reverse=True)
    return all_proteins


def profile_matrix(seq_dict):
    '''
    Function to return the profile matrix for a collection of sequences
    Args:
        seq_dict: A dictionary containing the sequences
    Return:
        The profile matrix with row labels: [A, C, G, T]
    '''
    pMatrix = {}
    for seq in seq_dict:
        for i in range(len(seq)):
            if(seq[i] == "A"):
                if(not 0 in pMatrix.keys()):
                    pMatrix[0] = [0] * len(seq)
                pMatrix[0][i] += 1
            elif(seq[i] == "C"):
                if(not 1 in pMatrix.keys()):
                    pMatrix[1] = [0] * len(seq)
                pMatrix[1][i] += 1
            elif(seq[i] == "G"):
                if(not 2 in pMatrix.keys()):
                    pMatrix[2] = [0] * len(seq)
                pMatrix[2][i] += 1
            elif(seq[i] == "T"):
                if(not 3 in pMatrix.keys()):
                    pMatrix[3] = [0] * len(seq)
                pMatrix[3][i] += 1
    return pMatrix


def consensus_sequence(profileMatrix):
    '''
    Function to return the consensus string from a dictionary of strings
    Args:
        seq_dict: dictionary of sequences
    Return:
        The consensus string
    '''
    maxBaseIndices = [None]*len(profileMatrix[0])
    for i in range(len(profileMatrix[0])):
        currMax = 0
        for j in range(4):
            if(max(currMax, profileMatrix[j][i]) > currMax):
                currMax = max(currMax, profileMatrix[j][i])
                maxBaseIndices[i] = str(j)
    return "".join(maxBaseIndices).replace("0", "A").replace("1", "C").replace("2", "G").replace("3", "T")


def motif_enumeration(dna_seqs, k, d):
    '''
    A function to return all (k,d)-motifs appearing in list of DNA sequences
    Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif
    if it appears in every string from Dna with at most d mismatches
    Args:
        dna_seqs: A list of valid DNA sequences
        k: length of motifs
        d: no. of max. mismatches
    Return:
        A list containing all (k,d)-kmer motifs
    '''
    kmer_motifs = []
    kmer_all = []
    for seq in dna_seqs:
        for i in range(len(seq) - k + 1):
            kmer_all.append(kmer_neighbours(seq[i:i + k], d))
    kmer_all = list(set([item for elem in kmer_all for item in elem]))
    for kmer in kmer_all:
        counter = 0
        for seq in dna_seqs:
            if len(motif_search_overlapping_approx(seq, kmer, d)) > 0:
                counter += 1
        if counter == len(dna_seqs):
            kmer_motifs.append(kmer)
    return set(kmer_motifs)


def score_of_kmer_in_dna_list(kmer, dna_list):
    '''
    Function to return the score d(Pattern, Dna)
    Given a k-mer Pattern and a set of strings Dna = {Dna1, … , Dnat},
    we define d(Pattern, Dna) as the sum of distances (min. hamming distance for all possible kmers in a dna string)
    between Pattern and all strings in Dna
    Args:
        kmer: a pattern
        dna_list: a list of dna strings
    Return:
        The score defined in function definition
    '''
    score = 0
    k = len(kmer)
    for dna in dna_list:
        d = len(kmer)
        for i in range(len(dna) - k + 1):
            if (hamming_distance(kmer, dna[i:i + k]) < d):
                d = hamming_distance(kmer, dna[i:i + k])
        score += d
    return score


def median_string(k, dna_list):
    '''
    Functio to return the median string
    median string is a k-mer Pattern that minimizes d(Pattern, Dna) over all k-mers Pattern
    Args:
        k: lenghth of pattern
        dna_list: a list of dna strings

    Return:
        A median string
    '''
    # set median to arbitrary kmer
    median = "A" * k

    # to get all possible kmer, choose an arbitrary kmer and find its d-neighbours where d = k
    all_possible_kmers = kmer_neighbours(median, k)

    # set distance to maximum possible score
    distance = k * len(dna_list)

    for kmer in all_possible_kmers:
        d = score_of_kmer_in_dna_list(kmer, dna_list)
        if d < distance:
            distance = d
            median = kmer

    return median


def profile_probability(kmer, profile_matrix):
    """
    A function to return the probability of a kmer based on a profile matrix
    Args:
        kmer: the kmer
        profile_matrix: the profile matrix
    Return:
        The probability of kmer computed w.r.t the profile matrix
    """
    matrix_row = {"A": 0, "C": 1, "G": 2, "T": 3}
    probability = 1
    for i in range(len(kmer)):
        probability = probability*profile_matrix[matrix_row[kmer[i]]][i]
    return probability


def profile_probable_kmer(seq, k, profile_matrix):
    '''
    Function to return the most probable kmer in seq based on profile_matrix
    Args:
        seq: the sequence to pick kmers from
        k: length of kmer
        profile_matrix: profile matrix
    Return:
        The profile-most-probable kmer
    '''
    probable_kmer = "A"*k
    probab = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmer_probab = profile_probability(kmer, profile_matrix)
        if kmer_probab > probab:
            probab = kmer_probab
            probable_kmer = kmer
    return probable_kmer


def longest_common_substring(dna_list):
    '''
    A function to return the shared motif (longest common substring among a list of dna sequences)
    Args:
        dna_list: a list of dna sequences
    Return:
        The longest common substring
    '''
    shared_motif = ""
    short_string = min(dna_list, key=len)
    max_k = len(short_string)
    for i in range(max_k, 0, -1):
        for j in range(max_k - i + 1):
            substring = short_string[j:i+j]
            count_substring = 0
            for dna in dna_list:
                if dna.count(substring) > 0:
                    count_substring += 1
            if count_substring == len(dna_list):
                shared_motif = substring
                break
        if len(shared_motif) > 0:
            break
    return shared_motif
