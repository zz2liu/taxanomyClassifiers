"""
9/27/2007 created from test_taxonomy_from_tipnames.py
"""
from util import RdpTaxonomy, Taxonomy
from py_util.unit_test import TestCase, main, set_trace
from py_util.tree_new import traverse, DndParser, lineage
from taxonomy.lineage_tree import build_genus_tree
from numpy import array
class TaxonomyTests(TestCase):
    def setUp(self):
        taxonomy_str = """((a, b)A, (c,d)B)root;"""
        self.obj = DndParser(taxonomy_str, Taxonomy)

    def test_getTaxon(self):
        obs = [(n.Name, n.getTaxon()) for n in traverse(self.obj)]
        self.assertEqual(obs, [('root', ['root']), ('A', ['root', 'A']),
            ('a', ['root', 'A', 'a']), ('b', ['root', 'A', 'b']),
            ('B', ['root', 'B']), ('c', ['root', 'B', 'c']), 
            ('d', ['root', 'B', 'd'])])

    def test_taxNodeFromTaxon(self):
        res = self.obj.taxNodeFromTaxon(['root', 'A', 'a'])
        self.assertEqual(res.Name, 'a')

        res = self.obj.taxNodeFromTaxon(['Root', 'a'])
        self.assertEqual(res.Name, 'A')
        set_trace()

class RdpTaxonomyTests(TestCase):
    def setUp(self):
        self.rdp3 = RdpTaxonomy.fromTaxLines(
                TEST_TAXA_STR.splitlines(), maxlevel=3) #maxlevel deprecated
        #self.p(self.rdp3._root)

        self.rdp5 = RdpTaxonomy.fromTaxLines(
                TEST_TAXA_STR.splitlines(), maxlevel=5)
        self.taxon_full = ['Bacteria', 'Aquificae', 'Aquificae',
                'Aquificales', 'Aquificaceae', 'Aquifex']
        self.tax_array_full = [0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0]
        self.tax_array_empty = [0] * len(self.tax_array_full)

    def test_indexTree(self):
        test = self.rdp3
        self.assertEqual(test._level_slices,
                [slice(0, 1, None), slice(1, 2, None), slice(2, 5, None),
                 slice(5, 8, None), slice(8, 11, None), slice(11, 15, None),
                 slice(15, None, None)])

    def test_taxArrayFromTaxon(self):
        test = self.rdp5
        full_length_and_more = self.taxon_full
        full_length = full_length_and_more[:5]
        half_length = full_length_and_more[:3]
        mismatch_halfway = full_length_and_more[:3] + ['Aquifii']
        mismatch_then_match = ['error'] + full_length_and_more[1:]
        match_mismatch_then_match = full_length_and_more[:3] + ['error']\
                + full_length_and_more[4:]
        with_subtax = full_length_and_more[:3] + ['sorangineae']\
                + full_length_and_more[3:]


        expect = self.tax_array_full
        empty = [0] * len(expect)
        expect_half = expect[:6] + empty[6:]
        self.assert_equal(test.taxArrayFromTaxon(full_length_and_more),
                expect)
        self.assert_equal(test.taxArrayFromTaxon(full_length), expect)
        self.assertEqual(test.taxArrayFromTaxon(half_length),
                expect_half)
        self.assertEqual(test.taxArrayFromTaxon(mismatch_halfway),
                expect_half)
        self.assertEqual(test.taxArrayFromTaxon(mismatch_then_match), empty)
        self.assertEqual(test.taxArrayFromTaxon(match_mismatch_then_match),
                expect_half)
        self.assertEqual(test.taxArrayFromTaxon(with_subtax), expect)


    def test_taxonFromTaxArray(self):
        test = self.rdp5
        #level_slices = [(1, 2), (2, 5), (5, 8), (8, 11), (11, None)]
        self.p(test._level_slices)
        rborder3 = 8
        full = self.tax_array_full
        empty = self.tax_array_empty
        half = full[:rborder3] + empty[rborder3:]
        half_then_diverse = full[:rborder3] + [1] + full[rborder3:-1]

        expect_full = map(str.lower, self.taxon_full[:5]) #maxlevel=5
        expect_empty = []
        expect_half = expect_full[:3]
        
        self.assertEqual(map(str.lower, test.taxonFromTaxArray(full)),
                expect_full)
        self.assertEqual(map(str.lower, test.taxonFromTaxArray(empty)),
                expect_empty)
        self.assertEqual(map(str.lower, test.taxonFromTaxArray(half)),
                expect_half)
        self.assertEqual(map(str.lower, test.taxonFromTaxArray(half_then_diverse)),
                expect_half)


    def test_from_raw_tree(self):
        raw_tree = build_genus_tree(rdp_tax_tablines)
        res = RdpTaxonomy.fromRawTaxTree(raw_tree, maxlevel='genus')
        node = [n for n in res._root.iterNodes()
                if n.Name == '_test genus'][0]
        self.assertEqual([n.Name.lower() for n in lineage(node)],
           ['root', 'archaea', '_test phylum', '', '', '', '_test genus'])
        self.assertEqual(res._level_slices, [slice(0, 1, None), slice(1, 2,
            None), slice(2, 4, None), slice(4, 6, None), slice(6, 8, None),
            slice(8, 10, None), slice(10, None, None)])

        taxon_test = ['Archaea', 'Crenarchaeota', 'Thermoprotei',
                'Sulfolobales', 'Sulfolobaceae', 'Sulfolobus',
                'Sulfolobus acidocaldarius']
        taxarray_test = [0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0]
        
        taxarray = res.taxArrayFromTaxon(taxon_test)
        self.assertEqual(taxarray, taxarray_test)
        taxon = res.taxonFromTaxArray(array(taxarray_test))
        self.assertEqual(map(str.lower, taxon),
                map(str.lower, taxon_test[:-1])) #tax tree to genus


##########
# test_data
TEST_TAXA_STR = """domain Bacteria  (0/99183/0)  (selected/total/search matches)
        phylum Aquificae (0/558/0)
            class Aquificae (0/558/0)
                order Aquificales (0/558/0)
                    family Aquificaceae (0/529/0)
                        genus Aquifex (0/6/0)
                        unclassified_Aquificaceae (0/3/0)
                    family Incertae sedis 2 (0/29/0)
                        genus Balnearium (0/3/0)
                        genus Desulfurobacterium (0/2/0)
            unclassified_Aquificae (0/0/0)
        phylum Thermotogae (0/85/0)
            class Thermotogae (0/85/0)
                order Thermotogales (0/85/0)
                    family Thermotogaceae (0/85/0)
                        genus Thermotoga (0/34/0)
        phylum Genera_incertae_sedis_BRC1 (0/15/0)
            genus BRC1 (0/15/0)
            unclassified_Genera_incertae_sedis_BRC1 (0/0/0)"""

rdp_tax_tablines = """\
0*Root*-1*0*no rank
1*Archaea*0*1*domain
2*Crenarchaeota*1*2*phylum
22*_test phylum*1*2*phylum
3*Thermoprotei*2*3*class
4*_test subclass*3*4*subclass
64*Sulfolobales*4*4*order
65*Sulfolobaceae*64*5*family
66*Sulfolobus*65*6*genus
72*Sulfolobus solfataricus*66*7*species
67*Sulfolobus acidocaldarius*66*7*species
85*Sulfurisphaera*65*6*genus
86*Sulfurisphaera ohwakuensis*85*7*species
83*Stygiolobus*65*6*genus
900*_test genus*22*3*genus
""".splitlines()

if __name__ == '__main__':
    main()
