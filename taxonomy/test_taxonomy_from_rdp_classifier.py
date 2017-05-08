from taxonomy_from_rdp_classifier import parse_tip_taxa, _get_tip_lineage
from py_util.unit_test import TestCase, main, set_trace, StringIO
from util import RdpTaxonomy

class Test__get_tip_lineage(TestCase):
    def test_basic(self):
        inputs = [
        'gamma66;  ; Root; 100%; Bacteria; 100%; Proteobacteria; 75%; Gammaproteobacteria; 66%',
        'gamma24;  ; Root; 100%; Bacteria; 88%; Proteobacteria; 39%; Gammaproteobacteria; 24%',
        'rootonly; ; Root; 100%',
        ]
        expects = [
        ('gamma66', ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria'], [100.0, 75.0, 66.0]),
        ('gamma24', ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria'], [88.0, 39.0, 24.0]),
        ('rootonly', [], []),
        ]
        for inp, exp in zip(inputs, expects):
            self.assertEqual(_get_tip_lineage(inp), exp)
        
class Test_parse_tip_taxa(TestCase):
    def setUp(self):
        self.rdp_tree = 'RDP_taxonomy.tre'
        pass

    def test_basic(self):
        rdp_taxonomy = RdpTaxonomy(file(self.rdp_tree))
        res = list( parse_tip_taxa(
                StringIO(rdp_txt_test), rdp_taxonomy))
        self.assertEqual(res, default_exp)
        
    def test_default(self):
        res = list(parse_tip_taxa(StringIO(rdp_txt_test)))
        self.assertEqual(res, default_exp)

    def test_threshold(self):
        res = list(parse_tip_taxa(StringIO(rdp_txt_test), threshold=0.8))
        self.assertEqual(res, [('SBZAe693_clip_257_357', ['bacteria']),
                ('SBYO474_clip_257_357', ['bacteria']),
                ('SBYT2076_clip_257_357', ['bacteria', 'bacteroidetes']),
                ('SBYB1939_clip_257_357', ['bacteria']),
                ('C13_G10_clip_257_357',
                ['bacteria', 'firmicutes', 'clostridia', 'clostridiales']),
                ('SBYN1953_clip_257_357', ['bacteria', 'proteobacteria']),
                ('SBXY1029_clip_257_357', ['bacteria']),
                ('C18_a04_clip_257_357', ['bacteria', 'bacteroidetes']),
                ('C23_j13_clip_257_357', ['bacteria'])])



 





rdp_txt_test = """\
Classifier: Naive Bayesian rRNA Classifier Version 2.0, July 2007
Taxonomical Hierarchy: The nomenclatural taxonomy of Garrity and Lilburn, release 6.0
Query File: all_from_clip_257_357_unique.fa
Submit Date: Tue Jul 24 01:07:11 EDT 2007
Confidence threshold: 80%
Symbol - after a sequence name indicates the results are obtained using reverse complement of that query sequence.
Lineage: Root(4109)
Details:
SBZAe693_clip_257_357;  ; Root; 100%; Bacteria; 100%; Proteobacteria; 75%; Gammaproteobacteria; 66%; Salinisphaerales; 12%; Salinisphaeraceae; 12%; Salinisphaera; 12%
SBYO474_clip_257_357;  ; Root; 100%; Bacteria; 88%; Proteobacteria; 39%; Gammaproteobacteria; 24%; Cardiobacteriales; 8%; Cardiobacteriaceae; 8%; Suttonella; 8%
SBYT2076_clip_257_357;  ; Root; 100%; Bacteria; 100%; Bacteroidetes; 80%; Flavobacteria; 49%; Flavobacteriales; 49%; Flavobacteriaceae; 48%; Vitellibacter; 6%
SBYB1939_clip_257_357;  ; Root; 100%; Bacteria; 100%; Actinobacteria; 12%; Actinobacteria; 12%; Rubrobacteridae; 12%; Rubrobacterales; 12%; Rubrobacterineae; 12%; Rubrobacteraceae; 12%; Solirubrobacter; 11%
C13_G10_clip_257_357;  ; Root; 100%; Bacteria; 100%; Firmicutes; 89%; Clostridia; 86%; Clostridiales; 86%; Peptostreptococcaceae; 33%; Anaerococcus; 32%
SBYN1953_clip_257_357;  ; Root; 100%; Bacteria; 100%; Proteobacteria; 91%; Gammaproteobacteria; 76%; Oceanospirillales; 50%; Alcanivoraceae; 43%; Alcanivorax; 43%
SBXY1029_clip_257_357;  ; Root; 100%; Bacteria; 100%; Actinobacteria; 32%; Actinobacteria; 32%; Coriobacteridae; 13%; Coriobacteriales; 13%; Coriobacterineae; 13%; Coriobacteriaceae; 13%; Eggerthella; 6%
C18_a04_clip_257_357;  ; Root; 100%; Bacteria; 99%; Bacteroidetes; 90%; Bacteroidetes; 38%; Bacteroidales; 38%; Rikenellaceae; 29%; Marinilabilia; 17%
C23_j13_clip_257_357;  ; Root; 100%; Bacteria; 100%; Firmicutes; 79%; Clostridia; 73%; Clostridiales; 72%; Lachnospiraceae; 59%; Butyrivibrio; 49%
"""

default_exp = [('SBZAe693_clip_257_357',
  ['bacteria',
   'proteobacteria',
   'gammaproteobacteria',
   'salinisphaerales',
   'salinisphaeraceae']),
 ('SBYO474_clip_257_357',
  ['bacteria',
   'proteobacteria',
   'gammaproteobacteria',
   'cardiobacteriales',
   'cardiobacteriaceae']),
 ('SBYT2076_clip_257_357',
  ['bacteria',
   'bacteroidetes',
   'flavobacteria',
   'flavobacteriales',
   'flavobacteriaceae']),
 ('SBYB1939_clip_257_357',
  ['bacteria',
   'actinobacteria',
   'actinobacteria',
   'rubrobacterales',
   'rubrobacteraceae']),
 ('C13_G10_clip_257_357',
  ['bacteria',
   'firmicutes',
   'clostridia',
   'clostridiales',
   'peptostreptococcaceae']),
 ('SBYN1953_clip_257_357',
  ['bacteria',
   'proteobacteria',
   'gammaproteobacteria',
   'oceanospirillales',
   'alcanivoraceae']),
 ('SBXY1029_clip_257_357',
  ['bacteria',
   'actinobacteria',
   'actinobacteria',
   'coriobacteriales',
   'coriobacteriaceae']),
 ('C18_a04_clip_257_357',
  ['bacteria',
   'bacteroidetes',
   'bacteroidetes',
   'bacteroidales',
   'rikenellaceae']),
 ('C23_j13_clip_257_357',
  ['bacteria',
   'firmicutes',
   'clostridia',
   'clostridiales',
   'lachnospiraceae'])]

if __name__ == '__main__':
    main()