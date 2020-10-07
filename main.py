from bio_structures import *
from bioinformatics_tools import *
from algorithmic_tools import *


seq = "GGGCGGCTCG"


# Print statements to validate functions


# # dna_seq_validate()
# convert sequence to DNA sequence in all uppercase
#  or if the sequence is not a valid, get 0
#DNA_seq = dna_seq_validate(seq)


# # dna_freq_counter()
# print("Frequency Counts: ", dna_freq_counter(DNA_seq))
# print(str(dna_freq_counter(DNA_seq)["A"]) + " " +ṇ
#       str(dna_freq_counter(DNA_seq)["C"]) + " " +
#       str(dna_freq_counter(DNA_seq)["G"]) + " " +
#       str(dna_freq_counter(DNA_seq)["T"]))


# # rna_transcription()
#print("The transcribed RNA is", rna_transcription(seq))
#print("The complementary strand is:", dna_reverse_compliment(dna_seq_validate(seq)))


# # printing fibonacci series
# for i in range(25):
#     print("Fibonacci Series element", i+1, ":", fibonacci(n=i, k=1))
# print(fibonacci(n=4, k=3))


# # fasta_to_dict()
# print(fasta_to_dict('sample_fasta.txt'))


# # gc_content()
# print(gc_content(seq))
#print(gc_content_subset("ATGCATGCATGCATGCATGC", 4))


# # max gc content from fasta file
# new_dict = fasta_to_dict('sample_fasta.txt')
# gc_dict = {}
# for key, value in new_dict.items():
#     gc_dict[key] = gc_content(value)
# maxKey = max(gc_dict, key=gc_dict.get)
# print(maxKey)
# print(gc_dict[maxKey]*100)


# # hamming_distance()
# seq1 = "GAGCCTACTAACGGGAT"
# seq2 = "CATCGTAATGACGGCCT"
# print("Hamming distance:", hamming_distance(seq1, seq2))


# # rna translation
# print(rna_translation(seq))


# # motif_search_overlapping()
# seq = "GATATATGCATATACTT"
# sebseq = "ATAT"
# print(motif_search_overlapping(seq, sebseq))


# # generate_reading_frames()
# print(generate_reading_frames(seq))


# # proteins_from_seq()
# print(proteins_from_steq(['I', 'M', 'T', 'H', 'M', 'T',
#                             'Q', 'G', 'N', 'V', 'A', 'Y', 'I', '_']))


# # proteins_from_seq_orfs()
# seq = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
# for p in proteins_from_seq_orfs(seq, 0, len(seq), True):
#     print(p)

# # proteins_from_seq_orfs() from fasta
# dna_seq = fasta_to_dict('sample_fasta.txt')
# seq = dna_seq['Rosalind_0768']
# for p in set(proteins_from_seq_orfs(seq, 0, len(seq), True)):
#     print(p)
