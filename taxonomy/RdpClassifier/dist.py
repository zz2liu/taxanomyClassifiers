from numpy import *
def idx_to_ij(nrows):
    lefts = arange(nrows)[:-1].cumsum()
    print lefts
    def f(idx):
        i = lefts.searchsorted(idx, side='right')
        left = lefts[i-1]
        j = idx-left
        right = left+i-1
        row = range(left,idx) +range(idx+1, right)
        row += list(lefts[j:]+(i-1))

        print i, j, row#, col
        return i, j, row#, col

    return f

get_ij = idx_to_ij(5)
inputs = range(10)
for idx in inputs:
    get_ij(idx)

