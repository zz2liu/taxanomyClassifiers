"""multimer.py

provide functions with multimer words.
"""
from subprocess import Popen, PIPE
from numpy import array, zeros, disp
from itertools import count
from pdb import set_trace

from py_util.str_ import maketrans_except
from py_util.dict_ import calc_freqs
from py_util.cbook import cross_comb
from py_util.numpy_ import Array
from py_util.file_ import File
from py_util.tree_new import DndParser
from py_util.seq_ import parse_fasta
from py_util.maths.distances import dice_dist_lines

# make a tranlation table, to standard upper dna sequences.
invalid_to_spaces = maketrans_except('atcguATCGU', ' ')
invalid_to_spaces = invalid_to_spaces.upper().replace('U', 'T')


def iter_words(seq, size):
    """yield each overlapping fix-size multimer from a seq.

    - seq: a string of dna/rna sequences.
    - size: size of each word.

    Note: word are in uppercase and U->T.
    Note: chars other than ATCGU are excluded before making the words.
    """
    seq = seq.translate(invalid_to_spaces)
    for part in seq.split():
        for i in range(len(part)-size+1):
            yield part[i:i+size]

def unique_words(seq, size, **kw):
    """return a set of unique overlapping multimer from a seq."""
    return set(iter_words(seq, size, **kw))

def word_freqs(seq, size, **kw):
    """return the freqs of overlapping multimer from a seq."""
    return calc_freqs(iter_words(seq, size, **kw))

class MultimerClassifier(object):
    """

    Note: should perform best when each seq has the same length.
    todo: invalid_to_seq tobe implemented
    """
    def __init__(self, word_size, alphabet, invalid_to_sep=True):
        self._word_size, self._alphabet, self._invalid_to_sep = (
                word_size, alphabet, invalid_to_sep)
        self._words = words = \
                [''.join(w) for w in cross_comb([alphabet]*word_size)]
        self._empty_row = zeros(len(words), int)
        self._word_idxs = dict((w, i) for i, w in enumerate(words))


    def rowFromSeq(self, seq):
        """return a freq row from a seq.

        Used for classifying.
        """
        word_idxs = self._word_idxs
        result = self._empty_row.copy()
        for word in iter_words(seq, self._word_size):
            result[word_idxs[word]] += 1
        return result

    def _rowFromDict(self, dic):
        """return 1darray from a freq-dict."""

    def _dictFromRow(self, row):
        """return a freq-dict from a 1darray."""

    def matFromSeqs(self, seqs):
        """get a seq (row) by word (col) freq matrix.

        Used for training.
        """
        result = map(self.rowFromSeq, seqs)
        return array(result)

#####
# extension, not orthogonal for now
def freq_mat_from_fasta(fasta, word_size, 
        fasta_parser=parse_fasta, verbose=True):
    """return a freq matrix and labels from fasta lines.
    """
    mc = MultimerClassifier(word_size, alphabet='ATCG')
    freq_matrix = []
    labels = []
    counter = count()
    for label, seq in fasta_parser(fasta):
        if verbose: disp('\r%s: %s' % (counter.next(), label), linefeed=False)
        labels.append(label)
        freq_matrix.append(mc.rowFromSeq(seq))
    return array(freq_matrix), (labels, mc._words) #rownames, colnames
    
def tree_from_fasta(fasta, word_size, fasta_parser=parse_fasta, byfreqs=True,
        stdout=PIPE, distout=None, verbose=True):
    """return a relaxed nj tree from fasta lines.

    - fasta: lines or filename
    - stdout=PIPE: pass to Popen, can be open(filename, 'w')
    - distout=TempFileName: a optional distance output filename.
    
    Label willbe tipname.
    Dice distance matrix will be used to build tree.
    """
    # get labels and freq_matrix
    freq_matrix, (labels, words) = freq_mat_from_fasta(
            fasta, word_size, fasta_parser, verbose=verbose)

    # run clearcut to make the rnj tree
    if verbose: disp('make dist matrix file...')
    dist_lines = dice_dist_lines(freq_matrix, labels,
            byval=byfreqs, verbose=verbose)
    dist_file = File.fromLines(dist_lines, name=distout)
    
    if verbose: disp('clearcut running...')
    p = Popen('clearcut --distance', stdin=dist_file, stdout=stdout, shell=True)
    ret = p.wait()
    if ret: #failed
        raise RuntimeError('cleartcut failed')
    return p
