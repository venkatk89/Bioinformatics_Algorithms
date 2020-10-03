from bioinformatics_tools import *
from mathematical_tools import *


seq = "AAAACCCGGT"

# convert sequence to DNA sequence in all uppercase
#  or if the sequence is not a valid, get 0
DNA_seq = dna_seq_validate(seq)


# Print statements to validate functions


# print("The Sequence is:",  DNA_seq)
# print("Frequency Counts: ", dna_freq_counter(DNA_seq))
# print(str(dna_freq_counter(DNA_seq)["A"]) + " " +
#       str(dna_freq_counter(DNA_seq)["C"]) + " " +
#       str(dna_freq_counter(DNA_seq)["G"]) + " " +
#       str(dna_freq_counter(DNA_seq)["T"])) 

#print("The transcribed RNA is", rna_transcription(seq))
#print("The complementary strand is:", dna_strand_compliment(dna_seq_validate(seq)))

# for i in range(10):
#     print("Fibonacci Series element", i, ":", fibonacci(n=i, k=1))

#print(fibonacci(n=5, k=3))

print(decaying_fibonacci(n=82, k = 1, decay_rate=20))