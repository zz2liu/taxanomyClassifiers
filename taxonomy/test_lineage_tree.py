from lineage_tree import *
from py_util.unit_test import TestCase, main, set_trace
#??? how to avoid hard coded here ???
bergy78_tax_tab_file = '/home/zongzhi/Projects/misc/zongzhi/taxonomy/data/bergy7.8.taxid'

class MiscTest(TestCase):
    def setUp(self):
        self.tree = build_genus_tree(rdp_tax_tablines)

    def test_pad_mid_tax(self):
        res = pad_mid_missing_tax(['a', '', 'b', '', '', 'c', ''])
        self.assertEqual(res,  ['a', 'a_', 'b', 'b_', 'b__', 'c', ''])

    def test_build_tree(self):
        res = self.tree
        self.assertEqual(str(res), #with species removed
                "(((((((Sulfolobus,Sulfurisphaera,Stygiolobus)Sulfolobaceae)Sulfolobales)'_test subclass')Thermoprotei)Crenarchaeota,('_test genus')'_test phylum')Archaea)Root;")

    def test_dpcofg(self):
        res = lineage_to_dpcofg(self.tree, lineage_test1[1:])
        self.assertEqual(res, ['Archaea', 'Crenarchaeota', 'Thermoprotei',
                'Sulfolobales', 'Sulfolobaceae', 'Sulfolobus'])

        res = lineage_to_dpcofg(self.tree, lineage_test2[1:])
        self.assertEqual(res, ['Archaea', '_test phylum', '',
                '', '', '_test genus'])
        #None, None, '_genera_incertae_sedis'])

class Test_make_getter_directly(TestCase):
    def setUp(self):
        self.dpcofg = make_getter(open(bergy78_tax_tab_file))

    def test_Finder_unknown_at_last(self):
        #unknown at the last
        res = self.dpcofg('Bacteria;Dehalococcoides;'
                'Dehalococcoides_genera_incertae_sedis' .split(';'))                                                    
        self.assertEqual(res, #['bacteria', 'dehalococcoides'])
            ['Bacteria', 'Dehalococcoides', '', '', '',
                'Dehalococcoides_genera_incertae_sedis'])

    def test_Finder_subtax_in_the_middle(self):
        #subtax in the middle
        res = self.dpcofg('Bacteria;Proteobacteria;Deltaproteobacteria;'
                'Myxococcales;Sorangineae;Polyangiaceae;Sorangium'.split(';'))
        self.assertEqual(res, ['Bacteria', 'Proteobacteria', 'Deltaproteobacteria',
                'Myxococcales', 'Polyangiaceae', 'Sorangium'])

    def test_Finder_unknown_in_the_middle(self):
        #unknown in the middle
        res = self.dpcofg('Bacteria;Proteobacteria;__Deltaproteobacteria;'
                'Myxococcales;Sorangineae;Polyangiaceae;Sorangium'.split(';'))
        self.assertEqual(res, #['bacteria', 'proteobacteria'])
                ['Bacteria', 'Proteobacteria', '', '', '', ''])

    def test_Finder_missing_in_the_middle(self):
        #missing ranks in the middle
        res = self.dpcofg('Bacteria;Cyanobacteria;Cyanobacteria;Family 1.1;'
                'Microcystis'.split(';'))
        self.assertEqual(res, ['Bacteria', 'Cyanobacteria', 'Cyanobacteria',
                '', 'Family 1.1', 'Microcystis'])

    #mvd to parse_greengenes.py
    def _test_Finder_subtax_missing_in_the_middle(self):
        #missing ranks from greengenes-RDP
        res = self.dpcofg('Bacteria; Cyanobacteria; Cyanobacteria; Subsection 3;'
                ' Family 3.1; Microcoleus; otu_374'.split('; '))
        self.assertEqual(res, ['Bacteria', 'Cyanobacteria', 'Cyanobacteria',
                '', 'Family 3.1', 'Microcoleus'])
 
########
# test data
lineage_test1 = ['Root', 'Archaea', 'Crenarchaeota', 'Thermoprotei', 
        '_test subclass', 'Sulfolobales', 'Sulfolobaceae', 'Sulfolobus',
        'Sulfolobus acidocaldarius']
lineage_test2 = ['Root', 'Archaea', '_test phylum', '_test genus']
#node_id, node_name, parent_id, depth, rank_name
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
