"""test_taxonomy_from_groupnames.py -

02/02/08: get rid of dep on debug
"""
from __future__ import division
from taxonomy_from_groupnames import *
from py_util.unit_test import TestCase, main

from StringIO import StringIO
from numpy import zeros, array, all
from operator import itemgetter

from py_util.table import Table
from util import calc_taxa_recovery
from _test_taxonomy_from_groupnames import (TEST_TAX_STR, TEST_PHYLO_STR,
    FULL_RDP_TAX)


class FuntionTests(TestCase):
    def setUp(self): 
        self.tax = TreeFromRdpLines(StringIO(TEST_TAX_STR))
        
    def test_TreeFromRdpLines(self):
        tax = self.tax
        self.assertTrue(isinstance(tax, TreeNode))
        self.assertEqual(len(list(tax.traverse())),
                20)

    def test_get_self_and_ancestor_names(self):
        tax = self.tax
        first_tip = list(tax.tips())[0]
        names = get_self_and_ancestor_names(first_tip)
        self.assertTrue(len(names) > 3)
        
    def test_pacename_to_rdpname(self):
        pace_rdp = {
            #return as is
            'Alphaproteobacteria': 'Alphaproteobacteria',
            #mapping using the dict
            'Gemmimonas': 'Gemmatimonas', #special diff
            'Cystobacterineae': 'Myxococcales', #subclasses and suborders
            #rescue the first word
            'Sulfobacillus.2': 'Sulfobacillus',
            'Acidobacteria-1':  'Acidobacteria',
            'Micromonospora et al.': 'Micromonospora',
            #not rescued b/c there is [,/]
            'Peptococcaceae/Acidaminococcaceae':
                    'Peptococcaceae/Acidaminococcaceae',
            'Erythro-, Porphyrobacter': 'Erythro-, Porphyrobacter',
            #not rescued as in a dict already
            'Deinococcus-Thermus':   'Deinococcus-Thermus',
            }

        inputs = pace_rdp.keys()
        expects = pace_rdp.values()
        results = map(pacename_to_rdpname, inputs)
        for obs, exp in zip(results, expects):
            self.assertEqual(obs, exp)


class TaxonomyFromGroupNames_Tests(TestCase):
    def setUp(self):
        self.tax_tree = TreeFromRdpLines(StringIO(FULL_RDP_TAX))
        self.phil_tree = DndParser(StringIO(TEST_PHYLO_STR))
        self.tax_phil = TaxonomyFromGroupnames(self.tax_tree, self.phil_tree)

    def test_init(self):
        tax_phil = self.tax_phil
        tax_nodes = tax_phil._tax_nodes
        self.assertEqual(len(tax_nodes),
                1584)

        dups = [(name, len(nodes)) for name, nodes in tax_nodes.items()]
        dups.sort(key=itemgetter(1))
        self.assertEqual(map(itemgetter(1), dups[-3:]),
                [2, 2, 3]) ##Nitrospira:phylum, class, genus
        
        phylo_tree = tax_phil._phylo_tree
        self.assertEqual(len(list(phylo_tree.tips())),
                70)
        node_names = [n.Name for n in phylo_tree.traverse()]
        self.assertTrue(None not in node_names)
        self.assertTrue(all(not name.startswith('edge.')
            for name in node_names))

    def test_find_tax_node(self):
        t = self.tax_phil
        inputs = ['Sulfobacillus',
                'Synergistes',
                "'Bacillus cereus et al.'", #' stripped
                'Bacilli',
                ]
        for input in inputs:
            self.p(t.find_tax_node(input) is not None)

    def test__node_to_taxon(self):
        pass

    def test_eachGroupnamesTaxon(self):
        t = self.tax_phil
        xls_group_taxa(t, 'group_taxa.tmp')

    def test_eachTipTaxon(self):
        t = self.tax_phil
        xls_tip_taxa(t, 'tst1', selector=lambda x:x.Name.startswith('L'))


def xls_group_taxa(tax_phylo, outname, **kw):
    outfile = file(outname, 'w')
    outfile.write('%s\t%s\t%s\n'
            % ('Arb_group_names_from_root', 'mapped_pace_name',
                'RDP taxon:'+str(TITLES)))
    for groupnames,  mapped_name, taxon in \
            tax_phylo.eachGroupnamesTaxon(strict=False, **kw):
        outfile.write('%s\t%s\t%s\n'
                % ('; '.join(groupnames),
                    mapped_name, '; '.join(taxon)))
    outfile.close()

def xls_tip_taxa(tax_phylo, outname, **kw):
    outfile = file(outname, 'w')
    outfile.write('%s\t%s\t%s\t%s\n'
            % ('tip_name', 'RDP taxon:'+str(TITLES),
            'mapped_pace_name', 'group names from root'))
    for tip_name, taxon, mapped_name, groupnames in \
            tax_phylo.eachTipTaxon(return_others=True, strict=False, **kw):
        outfile.write('%s\t%s\t%s\t%s\n'
                % (tip_name, '; '.join(taxon),
                    mapped_name, '; '.join(groupnames)))
    outfile.close()


class RunTests:#(TestCase):
    def setUp(self):  
        self.tax_tree = TreeFromRdpLines(StringIO(FULL_RDP_TAX))
    
    def test_run_taxa_rick(self):
        rick_phylo_tree = DndParser(file(
                '/home/zongzhi/ArbData/454A_rick_w_bacteria.ntree'))
        rick_tax_phylo =  TaxonomyFromGroupnames(
                self.tax_tree, rick_phylo_tree)

        outname = '/home/zongzhi/Desktop/group_taxa.xls'
        xls_group_taxa(rick_tax_phylo, outname)

        outname = '/home/zongzhi/Desktop/tip_taxa.xls'
        xls_tip_taxa(rick_tax_phylo, outname,
                selector=lambda x:x.Name.startswith('454A_'))

       
    def _test_run_compare_between_full_and_clip(self):
        """make a file compare mice full and clipped taxa side by side"""
        taxa = self.taxa
        taxa.calcGroupnameTaxa()

        full_taxs =taxa.calcSeqnameTaxa(
            file('/home/zongzhi/mus_full.ntree'))
        clipped_taxs =taxa.calcSeqnameTaxa(file('mnk_hmn_mus/data/clip_r357_100_NAST.ntree'))

        #clipped_unique_taxs =taxa.calcSeqnameTaxa(file('/home/zongzhi/mus_r357_unique.ntree'))
        #seqname_representatives = seqname_representatives_from_unique_fasta(
        #        file('mnk_hmn_mus/data/_unique_ruth_r357_100_NAST.fasta'))
        #clipped_taxs = dict((name, clipped_unique_taxs[' '.join(rpr.split('_'))])
        #        for name, rpr in seqname_representatives.iteritems())

        ##stat the ratios of same observes at each level
        taxas_full = [tax for name, tax in sorted(full_taxs.items())]
        taxas_clip = [tax for name, tax in sorted(clipped_taxs.items())]

        print taxa.Titles
        same_ratios = calc_taxa_recovery(taxas_full, taxas_clip, len(taxa.Titles))
        print 'new_scoring:', same_ratios
        #same_ratios = tax_comparison(taxas_full, taxas_clip, len(taxa.Titles))
        #print 'old_scoring:', same_ratios
        #new_scoring: [ 1.          0.93379791  0.92116069  0.92132332
        #0.85991379 0.93617021 0.        ]
        #old_scoring: [1.0, 0.93379790940766549, 0.91739367502726277,
        #0.91724328812621092, 0.85499999999999998, 0.93617021276595747, None]

        #[('domain', 1.0), ('phylum', 0.93379790940766549), ('class', 0.92116068984396382), ('order', 0.92132332499304981), ('family', 0.85991379310344829), ('genus', 0.93617021276595747), ('species', 0.0)]
        #phylum_idx = taxa.Titles.index('phylum')
        #outfile = file('mnk_hmn_mus/full_clip_taxa.xls', 'w')
        #outfile.write('Sequence_name\tFull_taxa\tClipped_taxa\n')
        #for (full_name, full_taxa), (clip_name, clip_taxa) in zip(
        #        sorted(full_taxs.items()), sorted(clipped_taxs.items())):
        #    if full_taxa[phylum_idx:phylum_idx+1] == clip_taxa[phylum_idx:phylum_idx+1]: continue
        #    outfile.write('%s\t%s\t%s\n' % (
        #        full_name.split('full')[0], '; '.join(full_taxa), '; '.join(clip_taxa)))
        #outfile.close()

    def _test_mnk454_taxs(self):
        """make a xls file with [mnk454 seqname, pace_tax_string, RDP taxa]"""
        taxa = self.taxa
        taxa.calcGroupnameTaxa()

        mnk454_rdptaxs = taxa.calcSeqnameTaxa(file('/home/zongzhi/r357_mnk.ntree'))
        mnk454_groupnames = taxa._seqname_taxnames

        outfile = file('mnk_hmn_mus/r357_mnk_taxs.xls', 'w')
        outfile.write('%s\t%s\t%s\n'
                % ('seqname', 'groupnames', '\t'.join(taxa.Titles)))
        for seqname in mnk454_rdptaxs.keys():
            rdptaxs = mnk454_rdptaxs[seqname]
            groupnames = mnk454_groupnames[seqname]
            outfile.write('%s\t%s\t%s\n'
                    % (seqname, '; '.join(groupnames), '\t'.join(rdptaxs)))
        outfile.close()

    
    def _test_compare_with_phil_taxonomy(self):
        result = _seqname_rdptaxs(file('454_monkey/master_annotation.tab'))
        print result[:10]




#################################
# TestData
if __name__ == '__main__':
    main()
