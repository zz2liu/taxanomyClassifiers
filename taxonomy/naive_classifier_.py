"""test_naive_classifier.py

08/22/07 create.

todo: convert prints to tests.
todo: test bootstrap?
"""
from naive_classifier import *
from py_util.unit_test import TestCase, main
from pdb import set_trace

class TestMisc(TestCase):
    def test_invalid_to_spaces(self):
        self.assertEqual('qatqqcguqU'.translate(invalid_to_spaces),
                ' AT  CGT T')
    def test_unique_words(self):
        #test upper, R->T, invalid, len(part)==size, len(part)<size.
        self.assertEqual(unique_words('atcgncuuraa', 3),
                set(['CTT', 'TCG', 'ATC']))
        #test empty
        self.assertEqual(unique_words('a', 3), set())
        self.assertEqual(unique_words('', 1), set())

    def test_bootstrap(self):
        self.p(list(bootstrap('12345', 3, 5)))
        #set_trace()


class TestClassifierI(TestCase):
    def setUp(self):
        self.obj_3 = NaiveClassifier(3)

class TestClassifier_train(TestClassifierI):
    Debug = False
    def setUp(self):
        super(TestClassifier_train, self).setUp()
    
    def _test__get_genus_seqs(self):
        obj = self.obj_3
        res = obj._get_genus_seqs(seq_genera, drop_single_seq_genera=True)
        self.p(res)
        self.assertEqualItems(res[0], [3,2])

        res = obj._get_genus_seqs(seq_genera, drop_single_seq_genera=False)
        self.p(res)
        self.assertEqualItems(res[0], [3,2,1])

    def _test__get_seq_counts(self):
        obj = self.obj_3
        genus_idxs_test = {'a': 0, 'b': 1, 'c': 2}
        counts, words = obj._get_seq_counts(seq_genera, genus_idxs_test)
        self.p(counts, words)
        self.assertEqual(len(counts), 9)
        self.assertEqual(counts[words['ACC']], [1,1,0])
        self.assertEqual(counts[words['CCC']], [3,2,1])

    def _test__get_word_priors(self):
        obj = self.obj_3
        seq_counts_test = array([[1, 1, 0],
               [1, 0, 0],
               [1, 0, 0],
               [1, 1, 0],
               [1, 0, 0],
               [1, 0, 0],
               [1, 0, 0],
               [1, 1, 0],
               [3, 2, 1]])
        genus_seqs_test = array([3,2,1])
        total_seqs_test = 6

        res = obj._get_word_posteriors(
                seq_counts_test, genus_seqs_test, total_seqs_test)
        self.p(res)
        self.assertFloatEqual(res, [
                [ 0.33928571,  0.45238095,  0.17857143],
                [ 0.30357143,  0.07142857,  0.10714286],
                [ 0.30357143,  0.07142857,  0.10714286],
                [ 0.33928571,  0.45238095,  0.17857143],
                [ 0.30357143,  0.07142857,  0.10714286],
                [ 0.30357143,  0.07142857,  0.10714286],
                [ 0.30357143,  0.07142857,  0.10714286],
                [ 0.33928571,  0.45238095,  0.17857143],
                [ 0.98214286,  0.97619048,  0.96428571]
                ])


    def test_train(self):
        obj = self.obj_3
        obj.train(iter(seq_genera), drop_single_seq_genera=False)
        self.pp(obj.__dict__)
        obj.train(seq_genera)
        self.pp(obj.__dict__)


class TestClassifier_classify:#(TestClassifierI):
    def setUp(self):
        super(TestClassifier_classify, self).setUp()

    def test__get_max_likelihood_genus(self):
        obj = self.obj_3
        words_tests = [
                ['ATA', 'CCC'],
                ['ACC', 'CCC'],
                ['CCC'],# should be 'c' but 'a'??
                ]

        for words_test in words_tests:
            self.pp(obj._get_max_likelihood_genus(words_test,
                    word_posteriors_test, word_idxs_test))

    def test__get_max_bootstrap_genus(self):
        obj = self.obj_3
        obj._word_size = 3
        obj._word_posteriors = word_posteriors_test
        obj._word_idxs = word_idxs_test

        seq_tests = [
               'ATGCCCCC',
               'ATACCCCC',
               'ATCCCCCC',
               'ATCCCCCC',
               'ACCCCCCC',
               'CCCCCCCC',
               ]
        for seq_test in seq_tests:
            self.pp(obj._get_max_bootstrap_genus(seq_test, 100))
        
    def test_classify(self):
        obj = self.obj_3
        obj._word_size = 3
        obj._word_posteriors = word_posteriors_test
        obj._word_idxs = word_idxs_test
        obj._genus_idxs = genus_idxs_test

        seq_tests = [
               'ATGCCCCC',
               'ATACCCCC',
               'ATCCCCCC',
               'ATCCCCCC',
               'ACCCCCCC',
               'CCCCCCCC',
               ]
        for seq_test in seq_tests:
            self.pp(obj.classify(seq_test, 100))


    def _test__classify(self):
        obj = self.obj_3
        obj._prob_distribution = prob_distribution_test
        self.pp(obj._classify(['ATC','CCC']))
        self.pp(obj._classify(['ACC','CCC']))
        self.pp(obj._classify(['CCC','CCC']))

#######
# test data
seq_genera = [ #(seq, genera)
        ('ATGCCCCC', 'a'),
        ('ATACCCCC', 'a'),
        ('ATCCCCCC', 'a'), #estimated to 'b'
        ('ATCCCCCC', 'b'),
        ('ACCCCCCC', 'b'),
        ('CCCCCCCC', 'c'), #estimated to 'a'
        ]

word_posteriors_test = array([
       [ 0.33928571,  0.45238095,  0.17857143],
       [ 0.30357143,  0.07142857,  0.10714286],
       [ 0.30357143,  0.07142857,  0.10714286],
       [ 0.33928571,  0.45238095,  0.17857143],
       [ 0.30357143,  0.07142857,  0.10714286],
       [ 0.30357143,  0.07142857,  0.10714286],
       [ 0.30357143,  0.07142857,  0.10714286],
       [ 0.33928571,  0.45238095,  0.17857143],
       [ 0.98214286,  0.97619048,  0.96428571],
       ])
word_idxs_test = {'ACC': 0,
      'ATA': 4,
      'ATC': 3,
      'ATG': 1,
      'CCC': 8,
      'GCC': 2,
      'TAC': 6,
      'TCC': 7,
      'TGC': 5}
genus_idxs_test = {'a': 0, 'b': 1, 'c': 2}


if __name__ == '__main__':
    main()
