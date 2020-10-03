from bioinformatics_tools import *
from algorithmic_tools import *


seq = "CCACCCTCGTGGTATGGCTAGGCATTCAGGAACCGGAGAACGCTTCAGACCAGCCCGGACTGGGAACCTGCGGGCAGTAGGTGGAAT"

# convert sequence to DNA sequence in all uppercase
#  or if the sequence is not a valid, get 0
#DNA_seq = dna_seq_validate(seq)

# Print statements to validate functions


# print("The Sequence is:",  DNA_seq)
# print("Frequency Counts: ", dna_freq_counter(DNA_seq))
# print(str(dna_freq_counter(DNA_seq)["A"]) + " " +ṇ
#       str(dna_freq_counter(DNA_seq)["C"]) + " " +
#       str(dna_freq_counter(DNA_seq)["G"]) + " " +
#       str(dna_freq_counter(DNA_seq)["T"]))

#print("The transcribed RNA is", rna_transcription(seq))
#print("The complementary strand is:", dna_strand_compliment(dna_seq_validate(seq)))

# # printing fibonacci series
# for i in range(10):
#     print("Fibonacci Series element", i, ":", fibonacci(n=i, k=1))

#print(fibonacci(n=5, k=3))
#print(decaying_fibonacci(n=6, k = 1, decay_rate=20))

# print(fasta_to_dict('sample_fasta.txt'))
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


# seq1 = "GAGCCTACTAACGGGAT"
# seq2 = "CATCGTAATGACGGCCT"
# print("Hamming distance:", hamming_distance(seq1, seq2))
