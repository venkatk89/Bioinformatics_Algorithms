# importing structures in bio_struct file
from bio_structures import *


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
