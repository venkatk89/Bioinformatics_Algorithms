from bioinformatics_tools import *

seq = "ATTCCGGA"

DNA_seq = dna_seq_validate(seq)

print("The Sequence is:",  DNA_seq)
print("Frequency Counts: ", dna_freq_counter(DNA_seq))


# print(str(dna_freq_counter(DNA_seq)["A"]) + " " +
#       str(dna_freq_counter(DNA_seq)["C"]) + " " +
#       str(dna_freq_counter(DNA_seq)["G"]) + " " +
#       str(dna_freq_counter(DNA_seq)["T"])) 


