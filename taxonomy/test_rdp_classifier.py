"""test_rdp_classifier.py
"""
from rdp_classifier import RdpClassifier, rdp_fasta_parser
from py_util.unit_test import TestCase, main, set_trace, hook_pdb
from _test_rdp_classifier import rdp_fasta_test, testQuerySeq_fasta
hook_pdb()


class TestMisc(TestCase):
    def test_rdp_fasta_parser(self):
        res = list(rdp_fasta_parser(rdp_fasta_test))
        self.p(res)

class TestRdpClassifier(TestCase):
    def setUp(self):
        self.rdp = RdpClassifier()

    def test_train(self):
        self.rdp.train(rdp_fasta_test)
        #self.pp(self.rdp.__dict__)

    def test_classify(self):
        rdp = self.rdp
        rdp.train(rdp_fasta_test, drop_single_seq_genera=False)
        result = list(rdp.classify(rdp_fasta_test))
        self.pp(result)

    def test_classify__testQuerySeq(self):
        rdp = self.rdp
        rdp.train(testQuerySeq_fasta, drop_single_seq_genera=False)
        result = list(rdp.classify(testQuerySeq_fasta))
        self.pp(result)


if __name__ == '__main__':
    main()



