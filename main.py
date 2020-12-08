from bio_structures import *
from bioinformatics_tools import *
from algorithmic_tools import *
import numpy as np


seq = "GGGCGGCTCG"


# Print statements to validate functions

# f = open("sample_fasta.txt", "r")
# for i in f:
#     dna_seqs.append(i.strip("\n"))

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
# seq1 = "GGGCCGTTGGT"
# seq2 = "GGACCGTTGAC"
# print("Hamming distance:", hamming_distance(seq1, seq2))


# # rna translation
# print(rna_translation(seq))


# # motif_search_overlapping()
# seq = "ACGACAGACGACAGACGACAGACGACAGGACACGACAGACGACAGGTTAACGACAGCCACGACAGCACGACAGACGACAGGTAACGACAGACGACAGACGACAGTACGACAGAGGACGACAGACGACAGACGACAGGACGACAGCACGACAGAAGCCACGACAGATAGTCACGACAGCACGACAGATACACGACAGCCACGAACGACAGAGAGACGACAGAGTTACGACAGAACGACAGGCACGACAGACGACAGACGACAGTACGACAGCTACGACAGCCAGGTACGACAGACGACAGTGATACGACAGCAGTGCACGACAGACGACAGACGACAGCCCAGCAAACGACAGACGACAGACGACAGGGGGCGGAGAACGACAGACGACAGCCCACACGACAGTCACGACAGCACGACAGGAGACGACAGAGCACGACAGAGAACTACGACAGTACGACAGAGGTTAAGAACGACAGACGACAGAACGACAGACCACGACAGGAAACGACAGTAATACGACAGGACGACAGGCTACGACAGTGACGACAGACGACAGCCGAGGTTTCACGACAGTGGGACGACAGAACGACAGACGACAGCACGACAGCACGACAGACGACAGGGGTCTTCGGACGACAGACGACAGTTATGAGGGACGACAGTAACGACAGAGATATACGACAGAACGACAGAAACGACAGACGACAGACGACAGGTCAGTATACGACAGGACGACAGAACGACAGCATTGAGAAACGACAGCCCTGATACGACAGACGACAGACGACAGGACGACAGGCACGACAGAACGACAGACGACAGTCTACGACAGATGGTACTACGACAGACGACAGAGACCACGACAGACGACAGTTCAACGACAGAACGACAGAGACGACAGACACGACAGACGACAGACACGACAGACGACAGGACACGACAGTCACGACAGTACACGACAGGGTCGACGACAGAAACGACAGACGACAGCTACGACAGATCCACGACAGGACGACAGTACCGTACGACAGAGACGACAGGGAGTACGACAGACGACAGACGACAGAGGACGACAGCCCACGACAGCTCGTCGGACGACAGACGACAGAACGACAGAATGACGACAGTCACGACAGGGGCACGACAGACGACAGACGACAGAGGTCACCGGAGTCTGCATTGGGGGCACCACGACAGCCACGACAGATACGACAGTACGACAGAACGACAGACGACAGCACGACAGGTTGGCGTTCACTACGACAGACGACAGACGACAGCAACACGACAGAATGACGACAGGACGACAGACGACAGGCCACGACAGACGACAGAACGACAGCGACGACAGTTTTACGACAGCGACGACAGACGACAGCTACGTCCATGTAACGACAGATCGACGACAGTGACGACAGACGACAGGCGGTGCTCACGACAGCTTACGACAGCCATACACGACAGACGACAGTTTCGACGACAGACGACAGACGACAGGACGACAGGTGGGGACGACAGCGTCAACGACAGACGACAGCTATGTACGACAGCTATGCTAACGACAGACGACAGGTTAGAACGACGACAGATACACGACAGTTACGACAGTGTTAAACGACAGACGACAGACGACAGGACGACAGACAATACGACAGACGACAGACGACAGCACAACGACAGCCACACGACAGTGACGACAGAAATACGACAGAACGACAGACGACAGGGTTGGGACGACAGGACGACAGAGACGACAGCAAGACGACAGACGACAGACGACAGGACGACAGTGACGACAGACGCACGACAGCACGACAGAAACGACAGGATATTACGACAGACGACAGTAACGACAGTGACGACAGACGACAGCAACGACAGAGTACGACAGTACGACAGGGTGTGTACGACGACAGTGTGACGGGAACGACAGTAACGACAGCACGACAGATCCACGGACGACAGGACGACAGAACCTCCACGACAGACGACAGAACGACAGCACGACAGCGAACGACAGACGACAGGCAACGACAGAACGACAGTACGACAGCCCACGACAGTCCACGACAGGGACGACAGACTACGACAGACGACAGGACGACAGACGACAGCTAACTTACGACAGACGACAGGTGTGATACACGACAGCACGCACGACAGACGACAGGAGACGACAGTCAAGGACGACAGCGACGACAGTTAACGACAGGACGACAGGGAACGACAGGGACGACAGTCACGCTACGACAGCACGACAGACGACAGTACGACAGGGACGACAGCGAACAAGTAGAGGATAATTGACGACAGAGGGACGACAGGGACGACAGCCACGACAGGAAGACGACAGTACGACAGCTACGACAGACGACAGATAGCGTGCACGACAGACGACAGTACGACAGACGACAGGACGACGACAGCACGACAGGACGACAGACGACAGACGACAGTACGACAGACGACAGACGACAGGCCACGACAGATGTACGACAGACGACAGGGACGACAGTTAGAGATCACGACAGGCACGACAGCGTGACGACAGAATACGACAGCCACGACAGACGACAGCACGACAGACGACAGACCACGACAGGGATTAGTCAAGGCAAAAACGACAGAACGACAGAAACGACAGACGACAGCGTGAAAACGACAGAGATTAAACGACAGAACGACAGGGCACGACAGACGACAGACGACAGCCACGACAGTCTTTACGACAGTACGACAGACGCCACGACAGCTAACGACAGCTGACGACAGACACGACAGCCACGACAGGTCGTATTACGACAGGCACACGACAGACGACAGCACGACAGACGACAGCCAACGACAGGACACGACAGACGACAGCTTTACGACAGAGCCCGCCGGTATAACCGACGACAGGACGACAGTCTACGACAGGACGACAGGCACGACAGGGACGACAGGGAACGACAGGGACGACAGGCTTCCACACCAATCGCAGCGACGACAGGGCACGACAGCACGACAGGGCTCAGAGCCTTTACGACAGAAACAAGCTACGACAGACGACAGAACGACAGATGACGACAGTGACGACAGTATACGACAGACGACAGCCGGACGACAGTTGACAGGTAACGACAGACGACAGTACGACAGGGTACGACAGCTTCTCACGACAGCCTTTACGACAGGGTACGACAGTGTACGACAGCCGTACGACAGATGCAGCACGACAGACGACAGGGACGACAGAACACGACAGACGACAGACGACAGAACGACAGCCAACTTACGACAGGACGACAGACGACAGGAACGACAGCCTTACAACACGACAGAGACGACAGACGACAGACGACAGTCACGACAGTCACGACAGTACGACAGACGACAGTACGACAGCCCTTGCGGCACGACAGACGACAGTGAAGCTAAACGACAGACGACAGCACGACAGACGACAGGACGACAGACGACAGAAACGACAGACGACAGCTGAACGACAGTACGACAGACACGACAGTGTACGACAGACGACAGTCCACTACGACAGTGTACGACAGCGCCAAACGACAGACGACAGGCGACGACAGTAACGACAGAACGACAGCACGACAGGGAACGACAGGCCACACGACAGGACCCGACGACGACAGTGACGACAGTACGACAGAATGTCCACGACAGTTAGTACACGACAGGGACTACCACGACAGGACGACAGCCGCCCGAGGGCCAACACCCCCAACGACAGATTGGGTACGACAGCGACGACAGACGACAGGAAACGACAGTCAACGACAGCACGACAGTACGACAGAAACGACAGACGACAGATCACGACAGCGGACGACAGCGGAACGACAGTCGACGACAGCCACGACAGACGACAGACGACAGCTCTCTACAAACGACAGTCTGTATACAGCACGACAGACGACAGACGACAGGACGACAGACGACAGCGACGACAGTTCCACAGAGACGACAGACGACAGAGATCACGACAGGACGTGTGACGACAGACGACAGGACTGATACGACAGGCGGCTGACACGACAGAACGACAGCCGACGACAGGCCACGACAGACACGACGACAGGACGACAGGAAACGACAGACGACAGCTCTACGGGGAGACGACAGCAATGCTGACGACAGACGACAGACGACAGACGACAGGTCAAACGACAGCTACGACAGACGACAGGAAAACGACAGCACGACAGTACGACAGCACGACAGAACGACAGACGACAGGTCAGTACGACAGTCACGACAGACGACAGCTCCAAACGACAGCATAAGCAACGACAGACGACAGTACACGACAGTACGACAGAACGACAGACGACAGTAGACGACAGACGACAGTACGACAGGTACGACAGAGACGACAGTGACGACAGTTGCCTTTACGACAGCATCACGACAGTTTCACGACAGCACTTTCTTTATGCTGCTTACGACAGACGACAGCGTACGACAGTGACGACAGAAGTCGCGCATGCCACGACAGGTCGCAAACGACAGGACGACAGTGTGCGAGCGACGACAGCGACGACAGGTGACGACAGCACGACAGGAAACGACAGTACGACAGTTAATAACGACAGTACGACAGAGACGGCCACGACAGACCACGACAGCGACGACAGACGACAGTCGTGAACGACAGTAACCGACGACAGCAACGCACGACAGCGACGACAGTAAAACGACAGACGACAGAACGACAGGGCACGACAGCGGGTACTACGACAGACGACAGAACGACAGATCGTGCACGACAGCCACAACGACAGACGACAGACGACGACAGGGAGACGACAGACGACAGTGCAGGATTACGACAGAAGGCGGAACGACAGACGACAGCACCCCTCCACGACAGGACGACAGTACGACAGACGACAGTAACTCACCGACGACAGAACGACAGACGACAGTCACGACAGTTGACGACAGAAAGACGACAGAGTAGACGACAGACGACAGAGACGACAGAACGACAGCACGACAGCTACGACAGACGACAGAACGACAGACTAAGCATTGTGATGTACGACAGACGACAGACGACAGGACGACAGATGTGATGACCCAGCCACGACAGAAGCTTACGACAGTTACGACAGAGCTACGACAGTACGACAGACGACAGACGACAGTCCTGGACGACAGACGACAGTACGACAGCTGGACGACAGTGACGACAGTCTAACGACAGCATTGACGACAGCGTCCTTACGACAGACGACAGCACGACAGACGACAGGCACACGACAGCTGCCCGACGACAGAGGGTCCTTCGATTACGACAGAAACGACAGACGACAGCACGACAGTTACGACAGGACGACAGACACACGACAGTTTATGATCGACACGACAGACGACAGAAACGACAGCACATACGACAGCACTGGCGTGTAACACGACAGCCGGACGACAGCACGACAGGGACGACAGTCTACGACAGTACGACAGAAACGACAGTAACGACAGGTTCAACGACAGACGACAGGACGACAGAACGACAGTATACGACAGGGTACGACAGACGACAGCGACGACACGACAGCCGGTAATTAACGACAGGCTACGACAGGATTACGACAGTAACGACAGCCTCAACGACAGACAACGACAGACGACAGTGAACGACAGACGACAGTGAGAGATTACGACAGACGACAGCTGAACGACAGCACGACAGACGGCATAACGACAGTGACGACAGGAGGATAACGACAGACACGACAGTTGGCTAGACGACAGACGACAGACGACAGGATGACGACAGGATTCCGACGACAGGGAATAGGGTTGACGACAGCATACGACAGACGACAGGACGACAGGCACGACAGGACGACAGAACGACAGACGACAGACGACAGTGGACGACAGACGACAGACAACGACAGTTGGACTACGACAGACGACAGACGACAGACGACAGGGGGCCACGACAGTACGACAGAGACGACAGTACGACAGGTTGTACGACAGTCTGGACGACAGTACGACAGACGACAGCTGAACGACAGACGACAGATCCCACGACAGTGTGTTACTTGCACGACAGTACGACAGTACGACAGATACGACAGCAACGACAGTACGACAGACGACAGGACGAGACGACAGACGACAGTACGACAGGGACGACAGGGCATTGGACGACAGCCACGACAGCGACGACAGTCCTGACGACAGGCCGTCACGACAGGTTTTGACGACAGGCTTCGAACGACAGACGACAGCAACTTGCAGACGACAGACGACAGGACGACAGCCGACTGACGACAGACGACAGGACGACAGACGACAGAAACGACAGAAACGACAGACGACAGCACGACAGACGACAGACGCACGACAGCCTCACGACAGTATACGACAGTACGACAGACGACAGTGAACGACAGGTACGACAGTATCAACGACAGAACAGCTGGAGCACGACAGTGAATAGCTACGACAGATGACGACAGACGACAGCTCAGATGCACGACAGAACGACAGACGACAGGCACGACAGGTGTACGACAGCACGACAGTCCGTACGACAGAGATCCTCTGGGTCGGTAGACGACAGGTGACGACAGCCACGGACGACAGTACATCGAACGACAGAACGACAGAGGCATAACGACAGCCTCCCAAACGACAGAACGACAGGACACGACAGCCACACGACAGGATACACGACAGCTACGACAGACGACAGCACGACAGACACGACAGTACGACAGGCCTGCTACATTACGACAGGACGACAGTGCCTACCGCCACGACAGACGACAGGACGACAGTGAACGACAGACGACAGCATCTTACGACAGACGACAGTTACGACAGTACGACAGACGACAGATACGACAGGAACGACAGGACGACAGACGACAGGACGACAGACGACAGACGACAGATACGACAGACGACAGACCTTCGCACGACAGAACGACAGCACGACAGCAAAAAGACGACAGACGACAGAACGACAGACGACAGGAACGACAGACGACAGCGTACGACAGCAATCACGACAGTACGACAGAACGACAGCAACGACAGACGACAGACGACAGACGACAGAGAACGACAGCAACGACAGAAAAACGACAGCTAACGACAGCACGACAGACGACAGAACGACAGGTGACGACAGTACGACAGGAGACGACAGGTACGACAGGTACGACAGTGACGACAGCCACGACAGAACGACAGCACGACAGGACGACAGACGACAGGACGACAGAGAAGACGACAGACGACAGGACGACAGACGACAGGCTCACGACAGACGACAGTGCAACGACAGCACGACAGACCAGGACGACAGACGACAGAACGACAGAACGACAGCAAAAACGACAGACGACAGGACGACAGGGGTACGACAGCTGACGACAGTACGACAGGCGAGAACGACAGTCTGTAACGACAGTTGTACGACAGATCGGACACGACAGACGACAGTCTACGACAGACGACAGTACGACAGCATACGACGACAGAACGACAGAACGACAGCTACGACAGGCTACGACAGGAGGCGCACGACAGGCGACGACAGAGGGAAACGACAGACGACAGTACGACAGGGACGACAGCTTTGTCCATCACGACAGTGGCCGCCACGACAGGACGACAGACACGACAGAACGACAGGACGACAGGGCGTGAAACGACAGCACGACAGTACGACAGCAATACGACAGACGACAGTCCTATAGGAAAGGACGACAGAACTCACGACAGTGCCGGACCGGCATCACGACAGTCACGACGACAGGGATCTACGACAGGACTGACGACAGACGACAGACGACAGTTCCAGTTTACGACAGCGACGACAGACGACAGCACGACAGAACGACAGGACGACAGACGACAGGACGACAGAAGGAGACGACAGTTGAATGTAATACGACAGCAACGACAGGACGACAGAACGACAGACGACAGACTACGACAGGACGACAGCACGACAGCACACACGACAGTAACGACAGTATGACGACAGAACGACAGTACGACAGGGACGACAGACGACAGGTAGACGACAGGAAATACGCCTTGCACGACAGAAAACGACAGTACGACAGACGACAGTACGACAGCACGACAGGACGACAGGGTACGACAGCCCACGACAGCGCCGCTCACGACAGTTGACGACAGAGCACGACAGACGACAGGCGGCAGCACGACAGGCTTGACGACAGTCGCACTGGCTACGACAGCACGACAGACGACAGTTTACGACAGACGACAGATTTTCTCATTACGACAGAATGTACGACAGAGACGACCCGAACGACAGACTCTAATTTCAACGACAGCACGACAGGACGACAGACGACAGACGACAGTCGGATACGACAGCTCACGACAGCTACGACAGACGACAGTCGTGTATCTACGACAGACCGACAGGCACGACAGACGACAGACGACAGTCACGACAGATACACGACAGAATACGACAGATACGACAGACGACAGAAACGACAGCATTATGTCCCGGAAGGACGACAGCACGACAGCACGACAGACGACAGACGACAGTTTATTTACGACAGGACCTTACACGACAGCAAAACGACAGACGACAGGATCGTGACGACAGCACGACAGAACGGTATTACGACAGACGACAGCCTCACGACAGACGACAGAGTACGACAGCACGACAGGACGACAGGTAACGACAGTTACACGACAGAACGACAGTCTGCAGAACGACAGACTGACGACAGACGACAGAACGACAGGACGACAGGACGACAGTCCTCCACGACAGAACGACAGACGACAGGAGCACGACAGAACGACAGCACGACAGACCTCGGACGACAGCCTACGACAGAGGTACGACAGGACGACAGTAACGACAGGCACGACAGACGACAGACGACAGTACATCACGACAGTGGTACGACAGAGGCTACGACAGCAAACGACAGACGACAGATACGACAGTGTAAGGATACGACAGGCACGACAGCACGACAGTAAACGACAGACGACAGACGACAGGGGGTGACGACAGAATATTCCTACGACAGTTAATCACGACAGGAAGCAACGGACGACAGCACGACAGCACGACAGACGACAGACGACAGTTTGCCTTACCAACGACAGTCGACGACAGGGTAACGACAGTGTTACGACAGACGACAGGACGACAGTCTAAACGACAGTCCACGACAGCGGTCACGACAGTACAACGACAGACGACAGCACGACAGACGACAGCACGACAGACGACAGACGACAGAACGACAGACGACAGATGATAAACCTTGGTGCAGACGACAGGCCCCAACGACAGCACGACAGACGACAGGACGACAGACGACAGACGACAGTCAGGGCTTACGACAGACGACAGTGGGACGACAGACGACAGCACGACAGGACGACAGTCACGACAGTGCGTTGACGACAGACGACAGACTGGACGCTTGACGACAGACGACAGCGAACGACAGACGACAGGACTACAACGAGATCTACGACAGTACGACAGCCTTGACGACAGTGCTGACGACAGAAACGACAGTTACGACAGTATGAAGAACGACAGGGACGACAGATACGACAGACGACAGACGACAGGACGACAGCATGTAGTGACACGACAGGCGTGGTACGACAGACGACAGACGACAGTACGACAGCCGGACACGACAGACGACAGACGACAGTTGACGACAGGCACGACAGACGACAGGAGACGACAGTGTCCGCACGACAGTGGACGACAG"
# sebseq = "ACGACAGAC"
# for i in motif_search_overlapping(seq, sebseq):
#     print(i-1, end=" ")


# # generate_reading_frames()
# print(generate_reading_frames(seq))


# # proteins_from_seq()
# print(proteins_from_steq(['I', 'M', 'T', 'H', 'M', 'T',
#                             'Q', 'G', 'N', 'V', 'A', 'Y', 'I', '_']))


# # proteins_from_seq_orfs()
# seq = 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG'
# for p in proteins_from_seq_orfs(seq, 0, len(seq), True):
#     print(p)


# # profile_matrix() and consensus_sequence()
# fDict = fasta_to_dict("sample_fasta.txt")
# seqList = fDict.values()
# print(consensus_sequence(profile_matrix(seqList)))
# for x in range(4):
#     print(Nucleotides[x], ": ", end="")
#     for j in profile_matrix(seqList)[x]:
#         print(j, end=" ")
#     print("\n")


# # most_frequent_kmer()
# text = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
# k = 4
# for i in most_frequent_kmer(text, k):
#     print(i, end=" ")


# # clump_finder()
# seq = 'CGGACTCGACAGATGTGAAGAAATGTGAAGACTGAGTGAAGAGAAGAGGAAACACGACACGACATTGCGACATAATGTACGAATGTAATGTGCCTATGGC'
# k = 5
# l = 75
# t = 4

# for i in clump_finder(seq, k, l, t):
#     print(i, end=" ")

# # min_skew_finder()
# seq = "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG"
# print(min_skew_finder(seq))


# # motif_search_overlapping_approx()
# seq = "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC"
# subseq = "ATTCTGGA"
# d = 3
# for i in motif_search_overlapping_approx(seq, subseq, d):
#     print(i, end=" ")


# # kmer_frequencies()
# seq = "ACGCGGCTCTGAAA"
# k = 2
# for i in kmer_frequencies(seq, k):
#     print(i, end=" ")


# # lexicographic_kmer
# print(lexicographic_kmer(7076, 11))


# # lexicographic_kmer_rank
# print(lexicographic_kmer_rank("TGTGCTGGAGAACTACCTATGCGGA"))

# # kmer_neighbours
# f = open("sample_fasta.txt", "w")
# kmer = "TAACTATCCT"
# d = 2
# for i in kmer_neighbours(kmer, d):
#     f.write(i)
#     f.write("\n")


# # most_frequent_kmer_approx
# seq = "AGTCAGTC"
# k = 4
# d = 2
# for i in most_frequent_kmer_approx(seq, k, d):
#     print(i, end=" ")


# # most_frequent_kmer_approx_reverse()
# seq = "ACGTTGCATGTCGCATGATGCATGAGAGCT"
# k = 4
# d = 1
# for i in most_frequent_kmer_approx_reverse(seq, k, d):
#     print(i, end=" ")


# # motif_enumeration
# k = 5
# d = 2
# dna_seqs = ['ATTTGGC', 'TGCCTTA', 'CGGTATC', 'GAAAATT']
# for i in motif_enumeration(dna_seqs, k, d):
#     print(i, end=" ")


# # median_string()
# dna_list = ["AAATTGACGCAT", "GACGACCACGTT",
#             "CGTCAGCGCCTG", "GCTGAGCACCGG", "AGTACGGGACAG"]
# k = 3
# print(median_string(k, dna_list))


# # profile_probability()
# kmer = "CCGAG"
# profile = [[0.2, 0.2, 0.3, 0.2, 0.3],
#            [0.4, 0.3, 0.1, 0.5, 0.1],
#            [0.3, 0.3, 0.5, 0.2, 0.4],
#            [0.1, 0.2, 0.1, 0.1, 0.2]]
# print(profile_probability(kmer, profile))


# # profile_probable_kmer()
# seq = "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT"
# k = 5
# profile = [[0.2, 0.2, 0.3, 0.2, 0.3],
#            [0.4, 0.3, 0.1, 0.5, 0.1],
#            [0.3, 0.3, 0.5, 0.2, 0.4],
#            [0.1, 0.2, 0.1, 0.1, 0.2]]
# print(profile_probable_kmer(seq, k, profile))
