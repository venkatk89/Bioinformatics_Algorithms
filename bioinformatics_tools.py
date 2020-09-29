
Nucleotides = ["A", "T", "C", "G"]


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
