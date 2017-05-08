from taxonomy_from_greengenes_classifier import *
from py_util.unit_test import TestCase, main, set_trace, StringIO

class Test_run_split_fasta(TestCase):
    def test_basic(self):
        split_fasta(file(split_fasta_test_filename), maxseqs=3)
        #ok


### test data
from py_util.seq_ import make_random_fasta
from py_util import TempFileName

split_fasta_test_filename = TempFileName(prefix='TempFileName_')
make_random_fasta(10, out_file=split_fasta_test_filename)



if __name__ == '__main__':
    main()
