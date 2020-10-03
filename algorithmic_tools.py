

def fibonacci(n, k = 1):
    '''
    Function to return the nth term of a fibonacci series
    Args:
        n: the value of n
        k: factor by which F(n-2) is multiplied with (no of rabbit pairs, one rabbits pair gives birth to)
    return:
        value: value of F(n)
    '''
    if n == 0:
        return 0
    if n == 1:
        return 1
    else: 
        return (fibonacci((n - 1), k) + k*fibonacci((n - 2), k))



def decaying_fibonacci(n, decay_rate, k = 1):
    '''
    Function to return the nth term of a fibonacci series
    Args:
        n: the value of n
        k: factor by which F(n-2) is multiplied with (no of rabbit pairs, one rabbits pair gives birth to)
        decay_rate: number of terms a rabbit pair would be alive
    return:
        value: value of F(n)
    '''
    if n <= 0:
        return 0
    if n == 1:
        return 1
    else: 
        return (decaying_fibonacci((n - 1), decay_rate, k) + k*decaying_fibonacci((n - 2), decay_rate, k) - decaying_fibonacci((n - decay_rate), decay_rate, k))