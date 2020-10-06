

# this algorithm is not computationally efficient
# def fibonacci(n, k = 1):
#     if n == 0:
#         return 0
#     if n == 1:
#         return 1
#     else:
#         return (fibonacci((n - 1), k) + k * fibonacci((n - 2), k))


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


# def decaying_fibonacci(n, decay_rate, k=1):
#     '''
#     Function to return the nth term of a fibonacci series
#     Args:
#         n: the value of n
#         k: factor by which F(n-2) is multiplied with (no of rabbit pairs, one rabbits pair gives birth to)
#         decay_rate: number of terms a rabbit pair would be alive
#     return:
#         value: value of F(n)
#     '''

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
