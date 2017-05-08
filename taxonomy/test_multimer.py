from multimer import *
from py_util.unit_test import TestCase, main, set_trace
import os


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

    def test_word_freqs(self):
        self.assertEqual(word_freqs('cuuut', 3),
                {'CTT': 1, 'TTT': 2})

class MultimerClassifierTests(TestCase):
    def setUp(self):
        self.mc = MultimerClassifier(2, 'ATC')
        
    def test__init(self):
        mc = self.mc
        self.assertEqual(mc._words,
                ['AA', 'AT', 'AC', 'TA', 'TT', 'TC', 'CA', 'CT', 'CC'])
        
    def test_rowFromSeq(self):
        res = self.mc.rowFromSeq('ATCCC')
        self.assertEqual(res, [0, 1, 0, 0, 0, 1, 0, 0, 2])

class MultimerClassifierTests_extension(TestCase):
    def test_freq_mat_from_fasta(self):
        res = freq_mat_from_fasta(fasta_test, 2)
        self.assertEqual(res[0].shape, (3, 16))
        
    def test_treeFromFasta(self):
        res = tree_from_fasta(fasta_test, 2)
        self.assertEqual(res.stdout.read().strip(),
                '((a:0.400000,b:0.000000),c:0.200000) ;')

        res = tree_from_fasta(fasta_test, 2, stdout=open('~tmp.tre', 'w'))
        self.assertEqual(open('~tmp.tre').read().strip(),
                '((a:0.400000,b:0.000000),c:0.200000) ;')
        os.remove('~tmp.tre')

#########
# test_data
fasta_test = """\
>a
ATTTTC
>b
ATTCCC
>c
CCCATT
""".splitlines()
if __name__ == '__main__':
    main()
