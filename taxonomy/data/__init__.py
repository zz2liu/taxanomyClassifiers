"""__init__.py

"""
from os import path

DIR = path.dirname(path.abspath(__file__))
BERGY51_TAX_TAB = path.join(DIR, 'bergy5.1.taxid')

BERGY78_TAX_TAB = path.join(DIR, 'bergy7.8.taxid')
BERGY78_DPCOFG = path.join(DIR, 'bergy7.8.dpcofg_list')
BERGY78_LINEAGES = path.join(DIR, 'bergy7.8_lineages.tab')
#use lineage_tree.make_getter to turn it into a (seqname: dpcofg) list.
#bergy5.1 taxonomy tree?
RDP_TREE_PATH = path.join(DIR, 'RDP_taxonomy.tre')
CORESET_LINEAGES = path.join(DIR, 'greengenes_coreset_rdp_tax.nds')
        #split by '; ' instead of ';'
        #62897 change from Aquificae to unclassified

TEST_PHYLO = path.join(DIR, '_test_from_groupnames.ntree')

bergy5_sub_taxs = [
        'subsection_1','subsection 1',
        'subsection_2','subsection 2',
        'subsection_3','subsection 3',
        'subsection_4','subsection 4',
        'Subsection_5','Subsection 5',
        'cystobacterineae',
        'sorangineae',
        'nannocystineae',
        'acidimicrobidae',
        'acidimicrobineae',
        'rubrobacteridae',
        'rubrobacterineae',
        'coriobacteridae',
        'coriobacterineae',
        'sphaerobacteridae',
        'sphaerobacterineae',
        'actinobacteridae',
        'actinomycineae',
        'micrococcineae',
        'corynebacterineae',
        'micromonosporineae',
        'propionibacterineae',
        'pseudonocardineae',
        'streptomycineae',
        'streptosporangineae',
        'frankineae',
        'glycomycineae']

bergy78_subclasses = ['caldilineae', 'actinobacteridae', 'acidimicrobidae',
       'coriobacteridae', 'sphaerobacteridae', 'rubrobacteridae']
bergy78_suborders = ['pseudonocardineae', 'streptomycineae',
        'propionibacterineae', 'streptosporangineae', 'micrococcineae', 'frankineae',
        'micromonosporineae', 'glycomycineae', 'corynebacterineae', 'actinomycineae',
        'acidimicrobineae', 'coriobacterineae', 'sphaerobacterineae',
        'rubrobacterineae', 'cystobacterineae', 'sorangineae', 'nannocystineae']
SUB_TAX_NAMES = dict.fromkeys(bergy5_sub_taxs + bergy78_subclasses +  bergy78_suborders)

