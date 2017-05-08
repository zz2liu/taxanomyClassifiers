"""test_eval_rdp_classifier.py
"""
from eval_rdp_classifier import *
from py_util.unit_test import TestCase, main, set_trace
from _test_rdp_classifier import testQuerySeq_fasta

class MiscTests(self):
    def test_parse_names_and_seqgenera(self):
        print parse_names_and_seq_genera(testQuerySeq_fasta)



if __name__ == '__main__':
    main()
