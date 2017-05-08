"""
9/27/2007 RdpTaxonomyTests mvd to test_util.py
"""
from __future__ import division
from pdb import set_trace
from numpy import (array, asarray, histogram, arange, ones, dtype, random, nan)
from pylab import hist, show
from profile import run as prun
from string import lower

from py_util.unit_test import TestCase, main
from py_util.tree_new import (DndParser, tips, nontips, traverse, lca,
        TreeNode, PhyloNode, search)
from taxonomy_from_tipnames import (PhyloTaxLCA, PhyloTaxLCA_LOO,
        TaxaTree, PhyloNode, TaxaNode,
        get_tipname_taxas, 
        tree_setKnownTaxas, tree_getUnknownTaxas, eval_taxa_result,
        _score_taxon, RdpTaxonomy, PhyloTaxTree, LeaveOneOutEvaluation)
from util import Taxonomy



class LcaTestI(TestCase):
    def setUp(self):
        taxonomy_str = """((a, b)A, (c,d)B)root;"""
        self.taxonomy = DndParser(taxonomy_str, Taxonomy)

        phylo_str = """(((1, 2)12, 3)123, 4)1234;"""
        self.phylo = DndParser(phylo_str, TreeNode)

        self.tipname_taxa = {'1': 'Aa', '2': 'Ab', '3': 'Bc', '4': 'Aa'}

        self.obj = PhyloTaxLCA(self.phylo, self.tipname_taxa, self.taxonomy, 
                _init_nodes=False)


class PhyloTaxLCATests(LcaTestI):
    def setUp(self):
        super(PhyloTaxLCATests, self).setUp()

    def test__init(self):
        obj = PhyloTaxLCA(self.phylo, self.tipname_taxa, self.taxonomy, 
                _init_nodes=False, taxon_with_root=False)

    def test__set_tip_tax_node(self): 
        obj = self.obj
        obj._set_tip_tax_node()
        tip_taxnodes = dict((tip.Name, tip.TaxNode.Name)
                for tip in tips(obj._root))
        self.assertEqual(tip_taxnodes, 
                {'1': 'a', '3': 'c', '2': 'b', '4': 'a'})

    def test__set_nontip_tax_node(self): 
        obj = self.obj
        obj._set_tip_tax_node()
        obj._set_nontip_tax_node()
        pnode_tnodes = dict((n.Name, n.TaxNode.Name)
                for n in nontips(obj._root))
        self.assertEqual(pnode_tnodes,
                {'1234': 'root', '12': 'A', '123': 'root'})

    def test_calc_tip_taxon_simple(self):
        obj = PhyloTaxLCA(self.phylo, self.tipname_taxa, self.taxonomy,
                taxon_with_root=False)

        obs = [(tip.Name, obj._calc_tip_taxon(tip))
                for tip in tips(obj._root)]
        self.assertEqual(obs, [('1', ['A']), ('2', ['A']),
                ('3', []), ('4', [])])

    def test_calc_tip_taxon_unknown(self):
        phylo_with_unknown_str = """\
                (
                    (
                        (
                            1,
                            x1,
                            (x2, 2)
                        )12,
                        3
                    )123,
                    4
                )1234;"""
        phylo = DndParser(phylo_with_unknown_str, TreeNode)
        obj = PhyloTaxLCA(phylo, self.tipname_taxa, self.taxonomy,
                taxon_with_root=False)
        nodes = dict((n.Name, n) for n in traverse(phylo))

        obs1 = obj._calc_tip_taxon(nodes['x1'])
        obs2 = obj._calc_tip_taxon(nodes['x2'])
        self.assertEqual(obs1, ['A'])
        self.assertEqual(obs2, ['A'])


class Tests_PhyloTaxLCA_LOO(LcaTestI):
    def setUp(self):
        super(Tests_PhyloTaxLCA_LOO, self).setUp()

    def test_getLineageRecoveries(self):
        obj = PhyloTaxLCA_LOO(self.phylo, self.tipname_taxa, self.taxonomy,
                taxon_with_root=False)
        obs = obj.getLineageRecoveries(levels=2, debug=True)
        #['root', None, None] ['root', 'a', 'a']
        #2 ['root', None, None] ['root', 'a', 'b']
        #3 ['root', 'A', None] ['root', 'b', 'c']
        #4 [None, None, None] ['root', 'a', 'a']
        self.assertEqual(obs,([[-1, -1], [-1, -1], [False, -1], [-1, -1]],
            ['1', '2', '3', '4']))

    
class Tests_PhyloTaxFitch_fitchback(TestCase):
    def setUp(self):
        taxonomy_str = """((a, b)A, (c,d)B)root;"""
        self.rdp = RdpTaxonomy(taxonomy_str)

        phylo_str = """(((1, 2)12, 3)123, 4)1234;"""
        self.phylo = DndParser(phylo_str, TreeNode)

        self.tipname_taxa = {'1': 'Aa', '2': 'Ab', '3': 'Bc', '4': 'Aa'}
        self.phylo_tax = PhyloTaxTree(self.phylo,
                self.tipname_taxa, self.rdp)

    def test__init(self):
        obj = PhyloTaxTree(self.phylo, self.tipname_taxa, self.rdp,
                fitch_back=True)
        obs = [(n.Name, obj._rdp.taxonFromTaxArray(n.TaxArray))
                for n in nontips(obj._root)]
        self.assertEqual(obs, [('1234', ['a', 'a']), 
                ('123', ['a', 'a']), ('12', ['a', 'a'])]) 

    def test__get_taxas_from_lines(self):
        inp = """a\tA; B;C\nb\t; BB; ;DD;\n""".splitlines()
        res = PhyloTaxTree._get_taxas_from_lines(inp)
        self.assertEqual(res,
                {'a': ['A', 'B', 'C'], 'b': ['', 'BB', '', 'DD', '']})


    def test_fitchback(self):
        obj = self.phylo_tax
        obs = [(n.Name, obj._rdp.taxonFromTaxArray(n.TaxArray))
                for n in nontips(obj._root)]
        self.assertEqual(obs, [('1234', ['a', 'a']),
                ('123', []), ('12', ['a'])])

        obj._adjust_nontip_tax_arrays(obj._root)
        obs = [(n.Name, obj._rdp.taxonFromTaxArray(n.TaxArray))
                for n in nontips(obj._root)]
        self.assertEqual(obs, [('1234', ['a', 'a']), 
                ('123', ['a', 'a']), ('12', ['a', 'a'])]) 

        #parent ambigous this time, no change
        obj = PhyloTaxTree(self.phylo.Children[0],
                self.tipname_taxa, self.rdp)
        obs = [(n.Name, obj._rdp.taxonFromTaxArray(n.TaxArray))
                for n in nontips(obj._root)]
        self.assertEqual(obs, [
                ('123', []), ('12', ['a'])])

        obj._adjust_nontip_tax_arrays(obj._root)
        obs = [(n.Name, obj._rdp.taxonFromTaxArray(n.TaxArray))
                for n in nontips(obj._root)]
        self.assertEqual(obs, [
                ('123', []), ('12', ['a'])])



class TestI(TestCase):
    def setUp(self):
        self.rdp = RdpTaxonomy.fromTaxLines(TEST_TAXA_STR.splitlines(), 5)
        self.phylo_tax = PhyloTaxTree(test_tree_str.splitlines(),
                test_nds_tabstr.splitlines(),
                self.rdp)
        self.rdp_tree_path = 'data/RDP_taxonomy.tre'
        a1, a2 = '88765', '88580'
        b1, b2 = '88790', '88640'
        c1 = '88635'
        self.aabbc_tree = """(
            (((xa:1, %s:1):1, %s:1):1, (xab:1, %s:1):1):1,
            (xb:1, (%s:1, (%s:1, xbc:1):1):1):1)root;""" % (
                a1, a2, b1, b2, c1)
        self.aabbc_root = DndParser(self.aabbc_tree.splitlines())
        self.aabbc = PhyloTaxTree(self.aabbc_root,
            test_nds_tabstr.splitlines(),
            file(self.rdp_tree_path))

class PhyloTaxTreeTests(TestI):
    """Tests for PhyloTaxTree.

    todo: get goot test data for _set_nontip_tax_arrays
    todo: modify the relevant tests where 
    """
    def setUp(self):
        super(PhyloTaxTreeTests, self).setUp()

    def test_init(self):
        test = self.phylo_tax
        self.p([n.Name for n in test._root.Children])
        self.p(test._taxas.items()[:5])
        self.p(test._rdp)

        #def test_setTipTaxArrays(self):
        #test = self.phylo_tax
        #test._init_tip_tax_arrays()
        tip_tax_arrays = [n.TaxArray for n in test._root.tips()]
        self.p(tip_tax_arrays[:3])
        #tip_taxons = [n._taxon for n in test._root.tips()]
        #self.p(tip_taxons[:3])

    def test_setNontipTaxArrays(self):
        test = self.phylo_tax
        #test._init_tip_tax_arrays()
        test._set_nontip_tax_arrays()
        tax_arrays = array([n.TaxArray for n in test._root.iterNodes()], int)
        print tax_arrays
        #set_trace()

    def test_calcTipTaxon(self):
        test = self.phylo_tax
        #test._init_tip_tax_arrays()
        test._set_nontip_tax_arrays()
        tip_taxons = [test._calc_tip_taxon(tip) for tip in test._root.tips()]
        self.p(tip_taxons[:10]) #the rdp taxonomy is too small for this test.

        test = self.aabbc
        test._set_nontip_tax_arrays()
        tip_taxons = dict((tip.Name, test._calc_tip_taxon(tip))
            for tip in test._root.tips())
        self.pp(tip_taxons)
        #set_trace()


class EvaluationTests(TestI):
    def setUp(self):
        super(EvaluationTests, self).setUp()
        self.ev = LeaveOneOutEvaluation(test_tree_str.splitlines(),
                test_nds_tabstr.splitlines(),
                file(self.rdp_tree_path))
        ## new test tree
        taxonomy_str = """((a, b)A, (c,d)B)root;"""
        self.rdp = RdpTaxonomy(taxonomy_str)

        phylo_str = """(((1, 2)12, 3)123, 4)1234;"""
        self.phylo = DndParser(phylo_str, TreeNode)

        self.tipname_taxa = {'1': 'Aa', '2': 'Ab', '3': 'Bc', '4': 'Aa'}

    def test_getLineageRecoveries(self):
        obj = LeaveOneOutEvaluation(self.phylo,
                self.tipname_taxa, self.rdp)
        obs = obj.getLineageRecoveries(levels=2, debug=False)
        #1 [None, None] ['a', 'a']
        #2 [None, None] ['a', 'b']
        #3 ['a', 'a'] ['b', 'c']
        #4 [None, None] ['a', 'a']
        self.assertEqual(obs, ([[-1, -1], [-1, -1], [False, False], [-1, -1]],
                ['1', '2', '3', '4']))
        obj = LeaveOneOutEvaluation(self.phylo,
                self.tipname_taxa, self.rdp, fitch_back=True)
        obs = obj.getLineageRecoveries(levels=2, debug=False)
        #1 ['a', None] ['a', 'a']
        #2 ['a', 'a'] ['a', 'b']
        #3 ['a', 'a'] ['b', 'c']
        #4 [None, None] ['a', 'a']
        self.assertEqual(obs, ([[True, -1], [True, False], [False, False], [-1,
            -1]], ['1', '2', '3', '4']))


    def test_score_taxon(self):
        f = LeaveOneOutEvaluation._score_taxon
        #exactly the same
        self.assertEqual(f('abc', 'abc'),
                1)
        #obs with more items
        self.assertFloatEqual(f('abcd', 'abc'),
                1)
        #obs with less items
        self.assertFloatEqual(f('ab', 'abc'),
                2/3)
        #obs empty
        self.assertFloatEqual(f('', 'abc'),
                0)
        #exp empty
        self.assertFloatEqual(f('abc', ''),
                nan)

        #obs wrong at beginning
        self.assertFloatEqual(f('bbcd', 'abc'),
                -1)
        #obs wrong at middle
        self.assertFloatEqual(f('abd', 'abcd'),
                -1/2)
        #obs wrong at end
        self.assertFloatEqual(f('abd', 'abc'),
                -1/3)

    def _test_calcScores(self):
        rdp = RdpTaxonomy(file(self.rdp_tree_path), maxlevel=5)
        tree = DndParser(file('/home/zongzhi/test_8k.ntree'))
        test = PhyloTaxTree(test_tree_str.splitlines(),
                test_nds_tabstr.splitlines(),
                rdp)
        scores, tipnames = test._calcScoresByLeftOneOut()
        n, bins = histogram(scores, bins=arange(-1,1,0.2))#, normed=True)
        self.p(n)

        scores, tipnames = test._calcScoresByLeftOneOut(exclude_leftout_tip=True)
        n, bins = histogram(scores, bins=arange(-1,1,0.2))#, normed=True)
        self.p(n)
        self.pp(sorted(zip(scores, tipnames))[:10])

    def _test_coreset_well_classified(self):
        rdp = RdpTaxonomy(file(self.rdp_tree_path), maxlevel=5)
        tree = DndParser(file('/home/zongzhi/coreset_well_classified.ntree'))
        nds_path = '/home/zongzhi/coreset_well_classified.nds'
        test = PhyloTaxTree(tree, file(nds_path), rdp)

        scores, tipnames = test._calcScoresByLeftOneOut()
        n, bins = histogram(scores, bins=arange(-1,1,0.2))#, normed=True)
        self.p(n)
        self.pp(sorted(zip(scores, tipnames))[:25])

        #scores, tipnames = test._calcScoresByLeftOneOut(exclude_leftout_tip=True)
        #n, bins = histogram(scores, bins=arange(-1,1,0.2))#, normed=True)
        #self.p(n)
        #self.pp(sorted(zip(scores, tipnames))[:10])

    def _test_external(self):
        rdp_tree_path = 'RDP_taxonomy.tre'
        rdp = RdpTaxonomy(file(rdp_tree_path))
        phylo_tax = PhyloTaxTree(test_tree_str.splitlines(),
                test_nds_tabstr.splitlines(), rdp)
        #phylo_tax._init_tip_tax_arrays()
        phylo_tax._set_nontip_tax_arrays()
        tip_taxons = array([phylo_tax._calc_tip_taxon(tip)
                for tip in phylo_tax._root.tips()], object)
        rand_idxs = random.permutation(len(tip_taxons))
        self.stop()
        self.p(tip_taxons[[rand_idxs[30]]])

    def _test_8k(self):
        rdp_tree_path = 'RDP_taxonomy.tre'
        rdp = RdpTaxonomy(file(rdp_tree_path))
        tree = DndParser(file('/home/zongzhi/test_8k.ntree'))
        phylo_tax = PhyloTaxTree(tree,
                file('/home/zongzhi/test_8k_.nds'), rdp)
        phylo_tax._init_tip_tax_arrays()
        phylo_tax._set_nontip_tax_arrays()
        tip_taxons = [phylo_tax._calc_tip_taxon(tip)
                for tip in phylo_tax._root.tips()]
        rand_idxs = random.permutation(len(tip_taxons))
        self.stop()
        self.p(tip_taxons[[rand_idxs[30]]])
        self.p(tip_taxons[:10])

    def _test_calcScores_8k(self):
        rdp_tree_path = 'RDP_taxonomy.tre'
        rdp = RdpTaxonomy(file(rdp_tree_path))
        tree = DndParser(file('/home/zongzhi/test_8k.ntree'))
        phylo_tax = PhyloTaxTree(tree,
                file('/home/zongzhi/test_8k_.nds'), rdp)
        scores, tipnames = phylo_tax._calcScoresByLeftOneOut()
        n, bins = histogram(scores, bins=arange(-1,1,0.1))#, normed=True)
        self.p(n)




class TaxaTreeTests:#(TestCase):
    def setUp(self):
        self.root = DndParser(test_tree_str, PhyloNode)
        self.taxas = TaxaTree._get_taxas_from_lines(test_nds_tabstr.splitlines())
        self.taxa_tree = TaxaTree(self.root, self.taxas)
        self.num_tips = len(list(self.root.iterTips()))

    def test_loadTree_from_node(self):
        test_from_node = TaxaTree()
        test_from_node.loadTree(self.root)
        self.assertEqual(test_from_node._root , self.root)

    def test_loadTree_from_lines(self):
        test_from_lines = TaxaTree()
        test_from_lines.loadTree(test_tree_str.splitlines())
        self.assertEqual(test_from_lines._root , self.root)

    def test_loadTaxas_from_dict(self):
        test_from_dict = TaxaTree()
        test_from_dict.loadTaxas(self.taxas)
        self.assertEqual(test_from_dict._taxas , self.taxas)

    def test_loadTaxas_from_lines(self):
        test_from_lines = TaxaTree()
        test_from_lines.loadTaxas(test_nds_tabstr.splitlines())
        self.assertEqual(test_from_lines._taxas , self.taxas)

    def _test_get_taxas_from_lines(self):
        self.assert_('NoTaxa' not in self.taxas)

    def test__initTreeTaxas(self):
        tree_with_taxas = self.taxa_tree
        tip_taxa_arrays = [t.TaxaArray 
                for t in tree_with_taxas._root.iterTips()]
        tip_taxas = [(t.Name, t.Taxa)
                for t in tree_with_taxas._root.iterTips()]

        #nontip_taxa_arrays = [n.TaxaArray 
        #        for n in tree_with_taxas._root.iterNontips()]

        self.assertEqual(len(tip_taxa_arrays), self.num_tips)
        self.assertEqual(array(tip_taxa_arrays).dtype, dtype(bool))
        self.assertEqualItems(array(tip_taxa_arrays).nonzero()[1],
                range(self.num_tips-1)) #one tip with NonTaxa
        #self.pp(asarray(tip_taxa_arrays[-5:], int),
        #        tip_taxas[-5:])
        #asarray(nontip_taxa_arrays[-5:], int))

    def test_calcNontipTaxaArrays(self):
        pass
    def test_calcNontipTaxa(self):
        pass

    def test_caleTipTaxa(self):
        taxa_tree = self.taxa_tree

        for tip in taxa_tree._root.tips()[:5]:
            self.p(taxa_tree.calcTipTaxa(tip, min_taxas=3), tip.Name, tip.Taxa)

            #    print taxa_tree.calcTipTaxa(tip, min_taxas=3), tip.Name, tip.Taxa
        #self.stop()


    
    def test_calcAllTipTaxas(self):
        tree_with_taxas = TaxaTree(self.root, self.taxas)

        #exclude tips with taxa
        test_12 = tree_with_taxas.calcAllTipTaxas(min_taxas=12)
        self.assertEqual(len(test_12), 1) # only one tip without taxa

        #include tips with taxa
        test_12 = tree_with_taxas.calcAllTipTaxas(exclude_taxa_tips=False,
                min_taxas=12)
        self.assertEqual(len(test_12), self.num_tips) # only one tip without taxa

class EvaluateTaxaTree:#(TestCase):
    def setUp(self):
        self.taxa_tree = TaxaTree(test_tree_str.splitlines(),
                test_nds_tabstr.splitlines())

    def test_evaluateByLeftOneOut_common(self):
        taxa_tree = self.taxa_tree
        scores, tipnames = taxa_tree._calcScoresByLeftOneOut(mode='common', 
                min_taxas=3)
        n, bins, patches = hist(scores, bins=arange(-1,1,0.1))#, normed=True)
        self.p(n)

        scores, tipnames = taxa_tree._calcScoresByLeftOneOut(mode='common',
                min_taxas=9)
        n, bins, patches = hist(scores, bins=arange(-1,1,0.1))#, normed=True)
        self.p(n)

    def test_evaluateByLeftOneOut_most(self):
        taxa_tree = self.taxa_tree

        scores, tipnames = taxa_tree._calcScoresByLeftOneOut(mode='most',
                min_taxas_most=9, threshold=0.9)
        n, bins, patches = hist(scores, bins=arange(-1,1,0.1))#, normed=True)
        self.p(n)

    def test_evaluateByLeftOneOut_both(self):
        taxa_tree = self.taxa_tree
        scores, tipnames = taxa_tree._calcScoresByLeftOneOut(mode='both',
                min_taxas=3, min_taxas_most=9, threshold=0.9)
        n, bins, patches = hist(scores, bins=arange(-1,1,0.1))#, normed=True)
        self.p(n)
        self.p([(score, tipname) for score, tipname in zip(scores, tipnames)
            if score < 0.1])
        #self.p(zip(scores, tipnames))
        #self.p(len(tipnames), len(set(tipnames)), len(scores))
    def test_evaluateByLeftOneOut_both(self):
        taxa_tree = self.taxa_tree
        f = lambda: taxa_tree._calcScoresByLeftOneOut(mode='both')
        self.stop()
        prun('f', 'tmp')

class FunctionTests:#(TestCase):
    def setUp(self):
        self.test_tree = DndParser(test_tree_str, TaxaNode)
        self.tipname_taxas = get_tipname_taxas(
                test_nds_tabstr.splitlines())
        self.tipname_taxas_part = get_tipname_taxas(
                test_nds_tabstr.splitlines()[:-10])
        pass

    def test_eval_taxa_result(self):
        #external test
        #tree = DndParser(file('/home/zongzhi/test_8k.ntree'),TaxaNode)
        #tipname_taxas = get_tipname_taxas(file('/home/zongzhi/test_8k_.nds'))
        def _eval(repeat=5, **kw):
            scores = []
            for i in range(repeat):
                result = eval_taxa_result(tree, tipname_taxas,
                        **kw)
                curr_scores = asarray(result, object)[:, 1]
                scores.extend(curr_scores.tolist())
            self.p(histogram(scores, bins=[-1, -0.5, 0, 0.2, 0.4, 0.6, 0.8]))
        #_eval(5, min_taxas=9)
        #_eval(5, min_taxas=18)
        #_eval(5, min_taxas=27)
        #_eval(5, min_taxas=27, threshold=0.9)
        #self.stop()

        tree = self.test_tree
        tipname_taxas = self.tipname_taxas
        scores = []
        for i in range(20):
            result = eval_taxa_result(self.test_tree, self.tipname_taxas)
            curr_scores = asarray(result, object)[:, 1]
            scores.extend(curr_scores.tolist())
        self.p(histogram(scores, bins=[-1, -0.5, 0, 0.2, 0.4, 0.6, 0.8]))

        _eval(20)
        _eval(20, min_taxas=12)
        _eval(20, min_taxas=15)
        _eval(20, min_taxas=20)
        _eval(20, min_taxas=20, threshold=0.9)

        #scores = []
        #for i in range(20):
        #    result = eval_taxa_result(self.test_tree, self.tipname_taxas,
        #            min_taxas=12)
        #    curr_scores = asarray(result, object)[:, 1]
        #    scores.extend(curr_scores.tolist())
        #self.p(histogram(scores, bins=[-1, -0.5, 0, 0.2, 0.4, 0.6, 0.8]))

        #scores = []
        #for i in range(20):
        #    result = eval_taxa_result(self.test_tree, self.tipname_taxas,
        #            min_taxas=15)
        #    curr_scores = asarray(result, object)[:, 1]
        #    scores.extend(curr_scores.tolist())
        #self.p(histogram(scores, bins=[-1, -0.5, 0, 0.2, 0.4, 0.6, 0.8]))

        #scores = []
        #for i in range(20):
        #    result = eval_taxa_result(self.test_tree, self.tipname_taxas,
        #            min_taxas=20)
        #    curr_scores = asarray(result, object)[:, 1]
        #    scores.extend(curr_scores.tolist())
        #self.p(histogram(scores, bins=[-1, -0.5, 0, 0.2, 0.4, 0.6, 0.8]))

        #scores = []
        #for i in range(20):
        #    result = eval_taxa_result(self.test_tree, self.tipname_taxas,
        #        min_taxas=20, threshold=0.9)
        #    curr_scores = asarray(result, object)[:, 1]
        #    scores.extend(curr_scores.tolist())
        #self.p(histogram(scores, bins=[-1, -0.5, 0, 0.2, 0.4, 0.6, 0.8]))
        self.stop()
    def test__get_taxa_score(self):
        self.assertEqual(_score_taxon('abc', 'abc'),
                1)
        self.assertFloatEqual(_score_taxon('ab', 'abc'),
                2/3)
        self.assertFloatEqual(_score_taxon('', 'abc'),
                0)

        self.assertFloatEqual(_score_taxon('b', 'abc'),
                -1)
        self.assertFloatEqual(_score_taxon('abd', 'abc'),
                -1/3)
        self.assertFloatEqual(_score_taxon('abcd', 'abc'),
                -1/4)

    def test_get_tipname_taxas(self):
        obs = get_tipname_taxas(test_nds_str.splitlines(), tab=False)
        self.p(obs.items()[:10])

        obs = get_tipname_taxas(test_nds_tabstr.splitlines(), tab=True)
        self.p(obs.items()[:10])

    def test_tree_setKnownTaxas(self):
        tree_setKnownTaxas(self.test_tree, self.tipname_taxas)
        self.pp([(t.Name, t.Taxa)
            for t in self.test_tree.tips()[:10]])

    def test_tree_getUnknownTaxas(self):
        tree_setKnownTaxas(self.test_tree, self.tipname_taxas_part)
        obs = tree_getUnknownTaxas(self.test_tree)
        self.pp(obs)


#################################
# Test Data
test_nds_tabstr = """\
NoTaxa	
8895	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Alteromonas
8890	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Alteromonas
8885	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Alteromonas
8880	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Alteromonas
8875	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Alteromonas
8870	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Alteromonadaceae; Glaciecola
8850	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Colwelliaceae; Colwellia
8840	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Colwelliaceae; Colwellia
8820	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Incertae sedis 7; Marinobacter
8815	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Incertae sedis 7; Marinobacter
8810	Bacteria; Proteobacteria; Gammaproteobacteria; Oceanospirillales; Hahellaceae; Hahella
8800	Bacteria; Proteobacteria; Gammaproteobacteria; Alteromonadales; Incertae sedis 7; Marinobacter
88995	Unclassified
88955	Bacteria; Bacteroidetes; Flavobacteria; Flavobacteriales; Flavobacteriaceae; Gillisia
88945	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Micrococcineae; Promicromonosporaceae; Isoptericola
88935	Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Hyphomicrobiaceae; Hyphomicrobium
88915	Bacteria; Proteobacteria; Alphaproteobacteria; unclassified_Alphaproteobacteria
88885	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Corynebacterineae; Nocardiaceae; Nocardia
88845	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Corynebacterineae; Mycobacteriaceae; Mycobacterium
88830	Bacteria; Firmicutes; Bacilli; Bacillales; Bacillaceae; Bacillus
88820	Bacteria
88815	Bacteria; Proteobacteria; Alphaproteobacteria; Caulobacterales; Caulobacteraceae; Caulobacter
88800	Bacteria; Firmicutes; Bacilli; Bacillales; Paenibacillaceae; Paenibacillus
88795	Bacteria
88790	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; unclassified_Comamonadaceae
88770	Bacteria; Proteobacteria; Gammaproteobacteria; Chromatiales; Chromatiaceae; unclassified_Chromatiaceae
88735	Bacteria; Proteobacteria; Alphaproteobacteria; Sphingomonadales; Sphingomonadaceae; Novosphingobium
88725	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Streptomycineae; Streptomycetaceae; Streptomyces
88710	Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Pseudomonadaceae; Pseudomonas
88700	Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Methylocystaceae; Methylocystis
88695	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Corynebacterineae; Gordoniaceae; Gordonia
88670	Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Hyphomicrobiaceae; Hyphomicrobium
88665	Bacteria; Proteobacteria; Alphaproteobacteria; Sphingomonadales; Sphingomonadaceae; Novosphingobium
88655	Bacteria; Proteobacteria; Alphaproteobacteria; Rhodobacterales; Rhodobacteraceae; Stappia
88640	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Incertae sedis 5; Aquabacterium
88635	Bacteria; Proteobacteria; Gammaproteobacteria; Pseudomonadales; Pseudomonadaceae; Pseudomonas
88630	Bacteria; Proteobacteria; Alphaproteobacteria; unclassified_Alphaproteobacteria
88595	Bacteria
88575	Bacteria; Proteobacteria; Betaproteobacteria; Rhodocyclales; Rhodocyclaceae; unclassified_Rhodocyclaceae
88570	Bacteria; Proteobacteria; Alphaproteobacteria; Sphingomonadales; Sphingomonadaceae; Sphingomonas
88560	Bacteria; Proteobacteria; Alphaproteobacteria; Sphingomonadales; Sphingomonadaceae; Sphingomonas
88510	Bacteria
88505	Bacteria; Proteobacteria; Alphaproteobacteria; Caulobacterales; Caulobacteraceae; Brevundimonas
88480	Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Bradyrhizobiaceae; unclassified_Bradyrhizobiaceae
88455	Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Arsenophonus
88445	Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Arsenophonus
88430	Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacteriales; Enterobacteriaceae; Arsenophonus
88425	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Hydrogenophaga
88400	Bacteria; Deinococcus-Thermus; Deinococci; Deinococcales; Deinococcaceae; Deinococcus
88360	Bacteria; Actinobacteria; Actinobacteria; Rubrobacteridae; Rubrobacterales; Rubrobacterineae; Rubrobacteraceae; unclassified_Rubrobacteraceae
88345	Bacteria; Proteobacteria; Alphaproteobacteria; Sphingomonadales; Sphingomonadaceae; Sphingomonas
88320	Bacteria; Cyanobacteria; Cyanobacteria; Subsection 4; Family 4.1; Aphanizomenon
88315	Bacteria; Spirochaetes; Spirochaetes; Spirochaetales; Serpulinaceae; Brachyspira
88300	Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Hyphomicrobiaceae; Hyphomicrobium
88295	Bacteria; Proteobacteria; Alphaproteobacteria; Rhizobiales; Hyphomicrobiaceae; Hyphomicrobium
88290	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Burkholderiaceae; Burkholderia
88260	Bacteria; Proteobacteria; Epsilonproteobacteria; Campylobacterales; Helicobacteraceae; Sulfurovum
88250	Bacteria; Proteobacteria; Alphaproteobacteria; Sphingomonadales; Sphingomonadaceae; Sphingomonas
88240	Bacteria; Firmicutes; Bacilli; Lactobacillales; Enterococcaceae; Enterococcus
88230	Bacteria
88205	Bacteria; Firmicutes; Bacilli; Lactobacillales; Streptococcaceae; Streptococcus
88200	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Hydrogenophaga
88180	Bacteria; Proteobacteria; Gammaproteobacteria; Oceanospirillales; unclassified_Oceanospirillales
88170	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; unclassified_Comamonadaceae
88130	Bacteria; Proteobacteria; Alphaproteobacteria; Caulobacterales; Caulobacteraceae; Brevundimonas
88115	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Comamonadaceae; Rhodoferax
88110	Bacteria; Firmicutes; Bacilli; Bacillales; Bacillaceae; Bacillus
88105	Bacteria; Proteobacteria; Betaproteobacteria; unclassified_Betaproteobacteria
88090	Bacteria; Proteobacteria; Betaproteobacteria; Burkholderiales; Burkholderiaceae; Pandoraea
88075	Bacteria; Firmicutes; Mollicutes; Incertae sedis 8; Erysipelotrichaceae; Erysipelothrix
88035	Bacteria; Bacteroidetes; Bacteroidetes; Bacteroidales; Prevotellaceae; Prevotella
88020	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Streptomycineae; Streptomycetaceae; Streptomyces
88155	Bacteria; Spirochaetes; Spirochaetes; Spirochaetales; Leptospiraceae; Leptospira
885	Unclassified
880	Unclassified
88805	Bacteria; Proteobacteria; Gammaproteobacteria; unclassified_Gammaproteobacteria
88765	Bacteria; Proteobacteria; Alphaproteobacteria; Caulobacterales; Caulobacteraceae; Brevundimonas
88745	Bacteria; Firmicutes; Bacilli; Lactobacillales; Lactobacillaceae; Lactobacillus
88715	Bacteria
88685	Bacteria; Actinobacteria; Actinobacteria; Actinobacteridae; Actinomycetales; Corynebacterineae; Nocardiaceae; Nocardia
88580	Bacteria; Proteobacteria; Alphaproteobacteria; Rhodospirillales; unclassified_Rhodospirillales
88350	Bacteria; Proteobacteria; Gammaproteobacteria; Xanthomonadales; Xanthomonadaceae; Thermomonas
88335	Bacteria; Actinobacteria; Actinobacteria; Rubrobacteridae; Rubrobacterales; Rubrobacterineae; Rubrobacteraceae; unclassified_Rubrobacteraceae"""

test_tree_str = """\
[The tree is based on tree_500_feb97. The partial sequences were added 
using the ARB_PARSIMONY tool. The topology of the initial tree was maintained.

[created as copy of 'tree_all'_]
(
  (
    (
      (
        (
          (
            (
              (
                (
                  (
                    (
                      (
                        (
                          (
                            (
                              (
                                (
                                  88455:0.00404,
                                  88445:0.00405
                                ):0.00000,
                                88430:0.00162
                              ):0.14375,
                              (
                                (
                                  (
                                    (
                                      (
                                        (
                                          8880:0.03583,
                                          (
                                            8875:0.03154,
                                            8885:0.01376
                                          ):0.00323
                                        ):0.00323,
                                        8895:0.02387
                                      ):0.02026,
                                      8870:0.07676
                                    ):0.01468,
                                    8890:0.02276
                                  ):0.03914,
                                  (
                                    8840:0.02273,
                                    8850:0.04702
                                  ):0.06165
                                ):0.03467,
                                NoTaxa:0.09974
                              ):0.00577
                            ):0.06235,
                            88180:0.17104
                          ):0.06955,
                          (
                            (
                              (
                                88635:0.07276,
                                88710:0.05990
                              ):0.07472,
                              (
                                8810:0.08864,
                                88805:0.05543
                              ):0.00327
                            ):0.00409,
                            (
                              (
                                8800:0.03325,
                                8820:0.02196
                              ):0.01212,
                              8815:0.01455
                            ):0.06976
                          ):0.00733
                        ):0.03755,
                        88770:0.08413
                      ):0.01140,
                      (
                        (
                          (
                            (
                              (
                                (
                                  (
                                    (
                                      88790:0.00161,
                                      88170:0.00250
                                    ):0.05177,
                                    88115:0.00325
                                  ):0.06879,
                                  (
                                    88200:0.00162,
                                    88425:0.00081
                                  ):0.04041
                                ):0.02111,
                                88640:0.12885
                              ):0.05271,
                              (
                                88090:0.04144,
                                88290:0.01291
                              ):0.04964
                            ):0.01547,
                            88575:0.08627
                          ):0.00407,
                          88105:0.14482
                        ):0.08194,
                        88350:0.15922
                      ):0.03854
                    ):0.01058,
                    (
                      (
                        (
                          (
                            (
                              88915:0.08480,
                              88630:0.09527
                            ):0.04494,
                            88655:0.03253
                          ):0.00574,
                          (
                            (
                              88480:0.19692,
                              88700:0.06877
                            ):0.04208,
                            (
                              (
                                (
                                  88295:0.00081,
                                  88300:0.00162
                                ):0.00325,
                                88935:0.00242
                              ):0.03592,
                              88670:0.06502
                            ):0.05139
                          ):0.00246
                        ):0.00568,
                        (
                          (
                            (
                              (
                                (
                                  88570:0.00000,
                                  88345:0.00081
                                ):0.00000,
                                88560:0.00575
                              ):0.04618,
                              88250:0.05923
                            ):0.00812,
                            (
                              88735:0.05592,
                              88665:0.00083
                            ):0.03643
                          ):0.12604,
                          88580:0.12838
                        ):0.01567
                      ):0.01053,
                      (
                        (
                          (
                            88505:0.00169,
                            88765:0.00085
                          ):0.05767,
                          88130:0.00575
                        ):0.07887,
                        88815:0.06660
                      ):0.12085
                    ):0.14702
                  ):0.05944,
                  88260:0.40548
                ):0.00989,
                (
                  88035:0.40393,
                  88955:0.19090
                ):0.23964
              ):0.00908,
              88510:0.23940
            ):0.00741,
            88400:0.22269
          ):0.01492,
          (
            (
              (
                (
                  88945:0.12910,
                  (
                    (
                      88845:0.06412,
                      88695:0.06744
                    ):0.02185,
                    (
                      (
                        88685:0.02679,
                        88885:0.00000
                      ):0.04201,
                      88995:0.00000
                    ):0.05569
                  ):0.06774
                ):0.00976,
                (
                  88020:0.07428,
                  88725:0.03393
                ):0.06414
              ):0.07648,
              (
                88335:0.09289,
                88360:0.09334
              ):0.04651
            ):0.04649,
            (
              (
                88320:0.22701,
                88230:0.11757
              ):0.05702,
              (
                (
                  88315:0.18859,
                  88155:0.16514
                ):0.03016,
                88595:0.22084
              ):0.01676
            ):0.01925
          ):0.00249
        ):0.01314,
        (
          88795:0.30326,
          88715:0.08246
        ):0.07127
      ):0.01310,
      (
        (
          (
            (
              (
                88745:0.19227,
                88240:0.03157
              ):0.01135,
              88205:0.18197
            ):0.05883,
            (
              88110:0.06489,
              88830:0.69240
            ):0.00081
          ):0.03587,
          88800:0.09130
        ):0.02487,
        88075:0.16478
      ):0.05193
    ):0.04835,
    88820:0.18342
  ):0.14987,
  (
    880:0.15812,
    885:0.10731
  ):0.12938
);"""
from test_util import TEST_TAXA_STR

if __name__ == '__main__':
    main()
