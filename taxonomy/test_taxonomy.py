from __future__ import division
#from __future__ import absolute_import
from StringIO import StringIO
from numpy import zeros, array

from py_util.debug import debug
from py_util.table import Table
from taxonomy_call import calc_taxa_recovery
from py_util.unit_test import TestCase, main
from taxonomy_ import *

def _seqname_rdptaxs(phil_tax_table):
    """get a dict of {seqname: rdp_taxs} from phil's master_annotation"""
    table = Table(phil_tax_table, header=True)
    seqnames_and_taxas = map(table.__getitem__,
            ['seqName', 'Phylum', 'Class', 'Order', 'Family', 'Genus'])
    result = {}
    for name_and_taxa in zip(*seqnames_and_taxas):
        name, taxa = name_and_taxa[0], name_and_taxa[1:]
        taxa = filter(None, taxa) # rid of '' from the end(assumption)
        result[name] = ['Bacteria',] + list(taxa)
    return result

class FuntionTests(TestCase):
    def test_misc(self):
        pass



def _firstone_firsttwo(line):
    fields = line.split()[:2]
    return fields[0], ' '.join(fields)

class TaxonomyTests(TestCase):
    def setUp(self):
        self.Empty = Taxonomy()

    def test_init(self):
        taxa = self.Empty
        self.assertEqual(taxa._taxa_tree, None)
        self.assertEqual(taxa._phylo_tree, None)

    def test_loadTaxaTree(self):
        taxa = self.Empty
        taxa.loadTaxaTree(StringIO(TEST_TAXA_STR))

        self.assertEqual(type(taxa._taxa_tree), TreeNode)
        print taxa.TaxaTree

        taxa.loadTaxaTree(StringIO(TEST_TAXA_STR), parser=_firstone_firsttwo)
        taxa.TaxaTree.pprint()
        #self.assertEqual(taxa._taxa_tree.NumDescendants, 10)

    def _test_loadPhyloTree(self):
        taxa = self.Empty
        taxa.loadPhyloTree(TEST_PHYLO_STR.splitlines()) #StringIO no works

        assert isinstance(taxa.PhyloTree, PhyloNode)


class RunTests:#(TestCase):
    def setUp(self):
        taxa = Taxonomy()
        taxa.loadTaxaTree(file('/home/zongzhi/ArbData/RDP_taxonomy.txt'))
        self.stop()
        taxa.loadPhyloTree(file('/home/zongzhi/phil16s_bacteria.ntree'))
        self.taxa = taxa
    
    def test_run_mapping(self):
        """make a file mapping pace taxa to RDP taxa"""
        taxa = self.taxa
        #a dict of {group_name: rdp_taxa }
        group_taxs, groupnames_taxs = taxa.calcGroupnameTaxa(verbose=True)
        outfile = file('mnk_hmn_mus/pace_RDP.xls', 'w')
        outfile.write('%s\t%s\n' % ('Pace taxa', 'RDP taxa'))
        for groupnames, taxs in sorted(groupnames_taxs.items()):
            outfile.write('%s\t%s\n' % ('; '.join(groupnames), '; '.join(taxs)))
        outfile.close()
       
    def test_run_compare_between_full_and_clip(self):
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
        debug.head(result)




#################################
# TestData
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

TEST_PHYLO_STR = ''


if __name__ == '__main__':
    main()
