from parse_rdp import *
from py_util.unit_test import TestCase, main, set_trace

class Test_each_lineage(TestCase):
    def test_basic(self):
        res = list(each_lineage(rdp_lines_test))
        #find all of them
        self.assertEqual(res, [
            ('SBXZ2962_clip_83_1326', ['Bacteria', 'Chlamydiae', 'Chlamydiae',
                'Chlamydiales', 'Parachlamydiaceae', 'Neochlamydia']),
            ('SBYS1055_clip_83_1326', ['Bacteria', 'Aquificae', 'Aquificae',
                'Aquificales', 'Aquificaceae', 'Hydrogenivirga']),
            ('SBYG1654_clip_83_1326', ['Bacteria', 'Cyanobacteria',
                'Cyanobacteria', 'Family 3.1', 'Arthrospira']),
            ('SBYX2890_clip_83_1326', ['Bacteria', 'Proteobacteria',
                'Deltaproteobacteria', 'Desulfobacterales',
                'Desulfobacteraceae', 'Desulfobacterium'])])
                     

    def test_min_percentage(self):
        res = list(each_lineage(rdp_lines_test, min_percentage=50))
        #find all of them
        self.assertEqual(res,   [
            ('SBXZ2962_clip_83_1326', ['Bacteria']),
            ('SBYS1055_clip_83_1326', ['Bacteria']),
            ('SBYG1654_clip_83_1326', ['Bacteria', 'Cyanobacteria',
                'Cyanobacteria']),
            ('SBYX2890_clip_83_1326', ['Bacteria', 'Proteobacteria',
                'Deltaproteobacteria', 'Desulfobacterales',
                'Desulfobacteraceae'])])

class Test_parse(TestCase):
    def test_basic(self):
        res = parse(rdp_lines_test, min_percentage=0)
        #find all of them
        self.assertEqual(res, [
            ('SBXZ2962_clip_83_1326', ['Bacteria', 'Chlamydiae', 'Chlamydiae',
                'Chlamydiales', 'Parachlamydiaceae', 'Neochlamydia']),
            ('SBYS1055_clip_83_1326', ['Bacteria', 'Aquificae', 'Aquificae',
                'Aquificales', 'Aquificaceae', 'Hydrogenivirga']),
            ('SBYG1654_clip_83_1326', ['Bacteria', 'Cyanobacteria',
                'Cyanobacteria', '', 'Family 3.1', 'Arthrospira']),
            ('SBYX2890_clip_83_1326', ['Bacteria', 'Proteobacteria',
                'Deltaproteobacteria', 'Desulfobacterales',
                'Desulfobacteraceae', 'Desulfobacterium'])])

    def test_threshold(self):
        res = parse(rdp_lines_test, min_percentage=50)
        self.assertEqual(res, 
            [('SBXZ2962_clip_83_1326', ['Bacteria', '', '', '', '', '']),
             ('SBYS1055_clip_83_1326', ['Bacteria', '', '', '', '', '']),
             ('SBYG1654_clip_83_1326', ['Bacteria', 'Cyanobacteria',
                 'Cyanobacteria', '', '', '']),
             ('SBYX2890_clip_83_1326', ['Bacteria', 'Proteobacteria',
                 'Deltaproteobacteria', 'Desulfobacterales',
                 'Desulfobacteraceae', ''])])
#mvd to lineage tree
#class Test_FromDpcofgTree:#(TestCase):

####
# test data
rdp_lines_test = """

Classifier: Naive Bayesian rRNA Classifier Version 2.0, July 2007
Taxonomical Hierarchy: Taxonomic Outline of the Bacteria and Archaea, release 7.8
Query File: gn__F343_100__unique.fa
Submit Date: Fri Dec 28 01:41:49 EST 2007
Confidence threshold: 80%
Symbol - after a sequence name indicates the results are obtained using reverse complement of that query sequence.
Lineage: Root(2509)
Details:
SBXZ2962_clip_83_1326;  ; Root; 100%; Bacteria; 100%; Chlamydiae; 26%; Chlamydiae; 26%; Chlamydiales; 26%; Parachlamydiaceae; 22%; Neochlamydia; 12%
SBYS1055_clip_83_1326;  ; Root; 100%; Bacteria; 100%; Aquificae; 37%; Aquificae; 37%; Aquificales; 37%; Aquificaceae; 37%; Hydrogenivirga; 35%
SBYG1654_clip_83_1326;  ; Root; 100%; Bacteria; 100%; Cyanobacteria; 50%; Cyanobacteria; 50%; Family 3.1; 24%; Arthrospira; 18%
SBYX2890_clip_83_1326;  ; Root; 100%; Bacteria; 100%; Proteobacteria; 93%; Deltaproteobacteria; 93%; Desulfobacterales; 84%; Desulfobacteraceae; 84%; Desulfobacterium; 8%
""".splitlines()
if __name__ == '__main__':
    main()
