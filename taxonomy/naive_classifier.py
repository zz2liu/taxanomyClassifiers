"""naive_classifer.py

Ref: Qiong Wang, et. al. 2007, AEM

08/22/07 created from the algorithm described.
done: unique_words ignore word with invalid bases
08/27/07 train default to be NOT drop_single_seq_genera.

todo: drop_single_seq_genera to seperate function.
todo: add lanemask and clip functions.
todo: _get_max_bootstrap always return a tuple, a solution to ties.
todo: some methods  to be functions or staticmethods.

Can be improved by assuming different generus prior prob.
Can be improved by taken word freqs (instead of unique words) into account?
"""
from __future__ import division
from collections import defaultdict
from operator import itemgetter
from itertools import count, tee

#from random import sample
from numpy import array, zeros, c_, concatenate

from py_util.dict_ import calc_freqs
from py_util.maths.math_ import sample
from taxonomy.multimer import unique_words, invalid_to_spaces

not_none = lambda x: x is not None


def bootstrap(items, choices, repeats):
    """yield each group of choices with replacement.

    - items: list like data.
    - choices: number of choices each time.
    - repeats: number of repeats.
    """
    for i in range(repeats):
        yield sample(items, choices, replace=True)


class NaiveClassifier(object):
    """holder o f train and classify.
    """
    def __init__(self, word_size=8, alphabet='ATCG'):
        """init."""
        self._word_size = word_size
        self._alphabet = list(alphabet)

    #####
    # train methods
    def train(self, seq_genera, drop_single_seq_genera=False):
        """generate a word ~ genera prob distribution matrix.

        - seq_genra: a seq of (seq, genus).
        set attributes the following attributes
        .word_posteriors: p(wi|G) in Formula 2. 
        .word_idxs: keep words in order.
        .genus_order: the unique genera as a list
        .genus_idxs: {genus: its idx in genus_order} for quick maping.

        Ref:
        Formula 1: p(wi) = (n(wi) + 0.5) / (N + 1)
        - p(wi): prob of observing the i th word in all sequences.
          (word_priors)
        - n(wi): numer of all sequences with the i th word. 
          (word_seqs)
        - N: numer of all the sequences in the training set. 
          (total_seqs)

        Formula 2: p(wi|G) = (m(wi) + p(wi)) / (M + 1)
        - p(wi|G): cond prob of observing the i'th word provided genus G.
          (word_posteriors)
        - M: number of sequences in genus G in the training set.
          (genus_seqs)
        - m(wi): number of sequences in G containing wi.

        """
        seq_genera, seq_genera_again = tee(seq_genera) #two iters
        #[#seqs in genus], {genus: idx}, #total_seqs.
        genus_seqs, genus_idxs, total_seqs = self._get_genus_seqs(
                seq_genera, drop_single_seq_genera)
        #a matrix of words ~ genera seqcounts and {word: idx}.
        seq_counts, word_idxs = self._get_seq_counts(
                seq_genera_again, genus_idxs)
        self._word_posteriors = self._get_word_posteriors(
                seq_counts, genus_seqs, total_seqs)
        self._genus_idxs, self._word_idxs = genus_idxs, word_idxs

    def _get_genus_seqs(self, seq_genera, drop_single_seq_genera):
        """return[#seqs in genus], {genus: idx}, total_seqs.
        """
        genera_with_repeats = [genus for seq, genus in seq_genera]
        genus_seqs = calc_freqs(genera_with_repeats)
        genera = array(genus_seqs.keys())
        counts = array(genus_seqs.values())

        if drop_single_seq_genera:
            mask = (counts > 1)
            genera = genera[mask]
            counts = counts[mask]

        genus_idxs = dict(zip(genera, count()))
        total_seqs = sum(counts)
        return counts, genus_idxs, total_seqs

    def _get_seq_counts(self, seq_genera, genus_idxs):
        """return a matrix of words ~ genera seqcounts and {word: idx}.

        Ref: m(wi) in Formula 1.
        """
        word_size = self._word_size
        #calc {word: [seqcounts for each genera]}
        result = defaultdict(lambda: zeros(len(genus_idxs), int))
        for seq, genus in seq_genera:
            try:
                genus = genus_idxs[genus]
            except KeyError, err: #genus dropped off
                continue
            for word in unique_words(seq, word_size):
                result[word][genus] += 1

        word_idxs = dict(zip(result.keys(), count()))
        seq_count_matrix = array(result.values())
        return seq_count_matrix, word_idxs

    def _get_word_posteriors(self, seq_counts, genus_seqs, total_seqs):
        """return a matrix of words ~ genera word_prob_provided_genus.

        use total_seqs as int, genus_seqs as int 1darray.

        Ref: Formula 1 and Formula 2.
        """
        # n(wi) as a col vector
        word_seqs = c_[seq_counts.sum(1)]
        # p(wi) as a col vector
        word_priors = (word_seqs + 0.5) / (total_seqs + 1)
        # p(wi|G)
        word_posteriors = (seq_counts + word_priors) / (genus_seqs + 1)
        return word_posteriors


    #######
    # classify methods
    def classify(self, seq, repeats=100):
        """return lineage from seq.

        genus_idxs to a common lineage.

        The joint prob of observing from genus G a sequence S, containing a set
        of words V(V <= W), was estimated as
        Formula 3: P(S|G) = prod( P(vi|G) )
        - P(vi|G): word posteria given by Formula 2.

        By Bayes' theorem, the prob that an unknown sequence, S, is a member of
        genus G is
        Formula 4: P(G|S) = P(S|G)P(G)/P(S)
        - P(S|G): joint prob given by formula 3
        - P(G) is the pior prob of S to be in G, assuming all genera are
          equally probable (equal priors), it becomes a constant. 
        - P(S) is the overall prob of observing S in any genera. Which can be
          ignored when classifing.

        So, classfy the seqence as a genus giving the highest joint prob.
        Formula 5: argmax( P(S|gi) )

        But, overlapping words in a query sequence is not independent.  So take
        only num_words/word_size words each bootstrap to calc the joint prob.
        Num_times to select a genus of 100 bootstrap trials was used as an
        estimate of confidence.  For higher-rank assignments, we sum for all
        genera under each taxon.
        """
        # calc genus_order
        genus_idxs = self._genus_idxs
        genera, idxs = map(array, (genus_idxs.keys(), genus_idxs.values()))
        genus_order = genera[idxs.argsort()]

        genus, bootstrap = self._get_max_bootstrap_genus(seq, repeats)
        return genus_order[genus], bootstrap

    def _get_max_bootstrap_genus(self, seq, repeats):
        """return the predicted genus idx and bootstrap value.

        - seq: a valid seqence as str.
        - repeats: repeats of bootstrap.
        use: unique_words, bootstrap, self
        """
        word_posteriors = self._word_posteriors
        word_idxs = self._word_idxs
        word_size = self._word_size

        all_words = list(unique_words(seq, word_size))
        print sorted(map(word_idxs.get, all_words))
        decisions = [] #genera idxs
        for words in bootstrap(all_words, len(seq)//word_size, repeats):
            decisions.append(self._get_max_likelihood_genus(words,
                    word_posteriors, word_idxs))
        freqs = calc_freqs(concatenate(decisions))
        sorted_freqs = sorted(freqs.items(), key=itemgetter(1))
        return sorted_freqs[-1] #what if a tie here?

    def _get_max_likelihood_genus(self, words,
            word_posteriors, word_idxs):
        """return the ML genus_idx from the words (features).

        """
        #Argmax prod( p(vi|G) )
        row_idxs = filter(not_none, map(word_idxs.get, words))
        likelihoods = word_posteriors[row_idxs].prod(0)
        # avoid .argmax() to solve tie problem.
        return (likelihoods == likelihoods.max()).nonzero()[0]


