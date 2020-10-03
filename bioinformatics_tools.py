
# Biological structures
Nucleotides = ["A", "T", "C", "G"]

Base_pairs = {"A":"T", "T":"A", "C":"G", "G":"C"}


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


def dna_strand_compliment(seq):
    '''
    Function to produce the reverse complement of a dna strand
    Args: The DNA sequence
    Return: The reverse compliment of the DNA strand. (BP is complimented, then string is reversed)
    '''
    temp_seq = ""
    for i in range(len(seq)):
        temp_seq += Base_pairs[seq[(len(seq) - 1 - i)]]
    return temp_seq

