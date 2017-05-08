from numpy import *

def idx_to_ij(nrows):
    """
    Eg, flat idxs of 5x5 dist matrix:
    [[],
    [0],
    [1,2],
    [3,4,5],
    [6,7,8,9]]
    """
    lefts = arange(nrows)[:-1].cumsum()
    def f(idx):
        i = lefts.searchsorted(idx, side='right')
        left = lefts[i-1] #left of curr row
        j = idx-left
        return i, j
    return f

get_ij = idx_to_ij(5)
inputs = range(10)
for idx in inputs:
    print get_ij(idx)

