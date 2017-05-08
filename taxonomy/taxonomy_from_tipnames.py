"""provide functions to get a dict of {seqname: RDP lineage}.

Classes provided:
    PhyloTaxLCA -> from last common ancestor of tax node.
    PhyloTaxFitch -> from fitch status of phylo node at each rank, with
        optional traverse back fitch assignments.

Terms used here:
    taxon: a list of classification names. eg. ['Bacteria', 'Proteo', ...]
    taxa: a list of taxons.
    tax: an item of taxon. eg. 'Proteo'
    taxs: a list of tax, but not forming a taxon. eg. ['Bacteria', 'Archea']
    Taxonomy: a hierarchy of taxons.
    Level: root is at level 0, root's children is at level 1 and so on.

> From greengenes.arb export nds {name > (acc, rdp tax)}
> insert alignments into the arb tree, export the tree.

07/23/2007 only use good ancestor to get lineage for a query tip.
12/04/2007 PhyloTaxLCA added.
02/02/08 mv calc_taxa_recovery to util.py


todo: add an evaluate method to show how good the tiptaxa have been recovered
    after searching the Taxonomy tree with the tipname_taxons.
todo: make the empty TaxArray as None to save space.
todo: reduce the dependency on the tree methods,  limited to:
    .Name, .Children, .Parent
    use functions for: traverse(), tips(), nontips(), ancestors()

todo: delete 'unclassified_*' items in rdp taxonomy?
todo: maxlevel, maxitems more consistent
todo: try BitSet instead of numpy int array for TaxArray
todo: move evaluate methods and tax recovery functions to another module.
todo: fromTaxLines move outof class, and merge with TreeFromLines in taxonomy_from_groupnames.
todo on PhyloTaxTree:
    Parrarel idea: only take a part of query seqs for one tree.

    Or only find the insert point to the barebone tree using parsimony insertion algorithm.
    If it works the classification algorithm will be really good.

    Another idea to improve the space complexity:
    build a rnj tree with db and query
    cut the query branches off, keeping the point nodes.
    
    Todo: RdpTaxonomy(Taxonomy)
"""
from __future__ import division
from pdb import set_trace
from numpy import (array, asarray, random, zeros, identity, logical_or,
    logical_and, all, any, nan, isnan)
from py_util.dict_ import count_freqs, resized #very light weight Freqs
from py_util.table import Table #simulation of data.frame in R
from py_util.tree_new import (PhyloNode, TreeNode, DndParser, lca, tips,
        nontips, search, ancestors, traverse)
from util import RdpTaxonomy

from taxonomy.lineage_recovery import get_lineage_recovery

or_reduce = logical_or.reduce
and_reduce = logical_and.reduce
permutation = random.permutation

def return_true(*a, **kw): return True
def notNone(x): return x is not None

class PhyloTaxLCA(object):
    """get taxon from the know tip-taxa by using last common ancester.

    A taxonomy tree and a phylo tree is needed.  Each tip of the phylo tree
    will be assigned to a node in the taxonomy tree first.  Then each nontip of
    the phylo tree will be assign to the last common ancestor of its children
    nodes. Each query tip inserted in the tree, should go back its ancestors to
    find a node with taxonomy assigned.
    """
    def __init__(self, phylo_tree, tipname_taxons, taxonomy, 
            _init_nodes=True, taxon_with_root=False):
        """
        - phylo_tree: the phylo tree built with db and query seqs together.
          as Newick lines or PhyloNode.
          using .tips(), .Name, .iterNontips(), .iterAncestors()
        - tipname_taxons: {tipname: taxon as str list}
        - taxon_with_root=False: if True, an input taxon should look like
          ['root', 'bacteria', ...]. output will be the similar
          
        - taxonomy: a Taxonomy Node. 
          using .taxNodeFromTaxon(), .getTaxon 
        _init_tax_arrays: for test only.
        """
        self._root = phylo_tree
        self._taxonomy = taxonomy
        self._taxas = tipname_taxons
        self._taxon_with_root = taxon_with_root
        if _init_nodes:
            self._set_tip_tax_node()
            self._set_nontip_tax_node()

    def _set_tip_tax_node(self):
        """assign a tax_node to each tip with taxon."""
        for tip in tips(self._root):
            taxon = self._taxas.get(tip.Name, None)
            if taxon:
                tip.TaxNode = self._taxonomy.taxNodeFromTaxon(
                        taxon, taxon_with_root=self._taxon_with_root)
            else:
                tip.TaxNode = None

    def _set_nontip_tax_node(self):
        for node in nontips(self._root, mode='postorder'):
            children_taxnodes = [c.TaxNode for c in node.Children
                    if c.TaxNode is not None]
            if children_taxnodes:
                node.TaxNode = lca(children_taxnodes)
            else:
                node.TaxNode = None

    def _calc_tip_taxon(self, tip):
        """return lineage for a tip or [].
        """
        for anc in ancestors(tip):
            children_known = [c for c in anc.Children
                    if c.TaxNode is not None]
            if len(children_known) > 1:
                return anc.TaxNode.getTaxon(self._taxon_with_root)
        else:
            return []

    def each_tip_taxon(self):
        """yield each (tipname, taxon calculated) for each unknown tip.
        """
        for tip in tips(self._root):
            if tip.TaxNode is None:
                yield tip.Name, self._calc_tip_taxon(tip)


class PhyloTaxLCA_LOO(PhyloTaxLCA):
    """leave one out evaluation of phylotax LCA method.
    """
    def getLineageRecoveries(self, levels=6, tip_filter=return_true,
            debug=False):
        """return [lineage recoveries] and corresponding [tipnames].

        - levels=6: start from root 
        """
        recoveries = []; tipnames = []
        for tip in filter(tip_filter, self._root.iterTips()):
            tipnames.append(tip.Name)
            #find and normalize the exp taxon for the curr tip
            exp_taxon = self._taxas.get(tip.Name, [])
            exp_taxon = [tax and tax.lower() or '' for tax in exp_taxon]
            exp_taxon = resized(exp_taxon, levels)

            #store and nullify TaxNode of the tip
            ori_TaxNode = tip.TaxNode
            tip.TaxNode = None
            #find the tip taxon after recalc the nontips
            self._set_nontip_tax_node()
            obs_taxon = self._calc_tip_taxon(tip) or []
            obs_taxon = resized(obs_taxon, levels)
            if debug: print tip.Name, obs_taxon, exp_taxon
            recoveries.append(get_lineage_recovery(obs_taxon, exp_taxon))
            #restore TaxNode of the tip
            tip.TaxNode = ori_TaxNode
        return recoveries, tipnames


class PhyloTaxTree(object):
    """get taxon for the tips by using Fitch status.
    
    PhyloTaxTree(tree, tipname_taxons, taxonomy)
    .calcAllTipTaxa(exclude_taxa_tips=True)

    Note: space complexity = k * 2*(m+n)  where
    k is number of end ranks(genura)
    m is number of db seqs
    n is number of query seqs
    """
    def __init__(self, phylo_tree, tipname_taxons, taxonomy, 
            _init_tax_arrays=True, fitch_back=False, **kw):
        """
        - phylo_tree: the phylo tree built with db and query seqs together.
          as Newick lines or PhyloNode.
        - tipname_taxons: {tipname: taxon as str list}
            or ['tipname\taxon split by ;']
        - taxonomy: a RdpTaxonomy or first arg of RdpTaxonomy()
        - fitch_back=False: preorder traversing, reduce the ambiguous status.
        **kw: pass to RdpTaxonomy(, maxlevel=Family)
        _init_tax_arrays: for test only.

        Note: the tipname_taxon is used to link the phylo_tree and taxonomy
        tree, and each taxon of it will be used to walk the taxonomy tree.
        """
        if not isinstance(phylo_tree, TreeNode):
            phylo_tree = DndParser(phylo_tree)
        self._root = phylo_tree

        if not isinstance(tipname_taxons, dict):
            tipname_taxons = self._get_taxas_from_lines(tipname_taxons)
        self._taxas = tipname_taxons

        if not isinstance(taxonomy, RdpTaxonomy):
            taxonomy = RdpTaxonomy(taxonomy, **kw)
        self._rdp = taxonomy
        self._fitch_back = fitch_back

        if _init_tax_arrays:
            self._init_tip_tax_arrays()
            self._set_nontip_tax_arrays()

    @staticmethod
    def _get_taxas_from_lines(name_taxa_lines):
        """return a dict of {tipname: taxa as a list}
        
        name_taxa_lines: each line is - tipname\ttaxon str sep by '; '
        """
        table = Table(name_taxa_lines, header=False)
        tipname_taxa = dict((name, map(str.strip, taxon.split(';')))
                for name, taxon in table
                if taxon) #not include in result if no valid taxa
        return tipname_taxa

    def _init_tip_tax_arrays(self):
        """
        """
        for tip in self._root.tips():
            taxon = self._taxas.get(tip.Name, [])
            tip.TaxArray = self._rdp.taxArrayFromTaxon(taxon).astype(bool)

    def _set_nontip_tax_arrays(self, start=None):
        if start is None: start = self._root
        self._set_nontip_tax_arrays_(start)
        if self._fitch_back:
            self._adjust_nontip_tax_arrays(start)

    def _set_nontip_tax_arrays_(self, start):
        """refresh .TaxArray for nontips from that of tips"""
        for node in start.iterNontips(mode='postorder'):
            node.TaxArray = self._get_tax_array_fitch(node)

    def _get_tax_array_fitch(self, node):
        """return the tax array for a nontip node by fitch method.
        
        algorithm: get children tax arrays; for each level, return intersection
        or union of it"""
        level_slices = self._rdp._level_slices
        children_tax_arrays = [c.TaxArray for c in node.Children]
        tax_array = and_reduce(children_tax_arrays)
        tax_array_union = or_reduce(children_tax_arrays)

        for level_slice in level_slices:
            if not any(tax_array[level_slice]):
                tax_array[level_slice] = \
                        tax_array_union[level_slice]
        return tax_array

    def _adjust_nontip_tax_arrays(self, start):
        """Tranversing preorder, if the state of the node is ambiquous, assign
        it to the parent if it is unambiguous.
        
        - start=self._root: the node to start adjusting.
        """
        each_node = traverse(start, 'preorder')
        root = each_node.next()
        for node in each_node:
            self._adjust_states(node)

    def _adjust_states(self, node):
        """if the state of the node is ambiquous, assign it to the
        parent if it is unambiguous."""
        level_slices = self._rdp._level_slices
        tax_array, parent_array = node.TaxArray, node.Parent.TaxArray

        for level in level_slices:
            if (tax_array[level].sum() > 1  #node ambiguous
                    and parent_array[level].sum() == 1): #parent unambiguous
                tax_array[level] = parent_array[level]
                
    def _calc_tip_taxon(self, tip):
        """return lineage for a tip or [].

        from its nearest good ancestor (derive from >= 2 known chidren).
        """
        for anc in tip.iterAncestors():
            children_known = [c for c in anc.Children
                    if any(c.TaxArray)]
            if len(children_known) > 1:
                return self._rdp.taxonFromTaxArray(anc.TaxArray)
        else:
            return []

    def each_tip_taxon(self):
        """yield each (tipname, taxon calculated) for each unknown tip.
        """
        for tip in tips(self._root):
            if not tip.TaxArray.any():
                yield tip.Name, self._calc_tip_taxon(tip)

    #deprecated, use each_tip_taxon instead
    def calcAllTipTaxa(self, exclude_taxa_tips=True,
            tiptaxon_getter=None):
        """return a dict of {tipname: taxon-predict}

        - exclude_taxa_tips=: if True, do not calc taxon for tips with taxon
        - tiptaxon_getter=self._calc_tip_taxon: a method to get taxon for a tip.
        """
        if not tiptaxon_getter:
            tiptaxon_getter = self._calc_tip_taxon

        result = {}
        for tip in self._root.iterTips():
            if exclude_taxa_tips and tip.TaxArray.any(): continue
            result[tip.Name] = tiptaxon_getter(tip)
        return result






class LeaveOneOutEvaluation(PhyloTaxTree):
    def _calc_tip_taxon_without_self(self, tip, debug=False):
        """return the taxon predicted for the tip.

        #zero tip.TaxArray, return,  then restore it.
        """
        ori_array = tip.TaxArray.copy()

        #nullify tip and recalc ancs
        tip.TaxArray[:] = 0
        self._set_nontip_tax_arrays()
        result = self._calc_tip_taxon(tip)

        #restore .TaxArray of self and ancs
        tip.TaxArray = ori_array
        return result

    def getLineageRecoveries(self, levels=6, tip_filter=return_true,
            debug=False):
        """return [lineage recoveries] and corresponding [tipnames].

        - levels=6: count from root if ._taxon_with_root, or count from 
          Domain.
        """
        recoveries = []; tipnames = []
        for tip in filter(tip_filter, self._root.iterTips()):
            tipnames.append(tip.Name)
            exp_taxon = self._taxas.get(tip.Name, [])
            exp_taxon = [tax and tax.lower() or '' for tax in exp_taxon]
            exp_taxon = resized(exp_taxon, levels)

            obs_taxon = self._calc_tip_taxon_without_self_fast(tip) or []
            obs_taxon = resized(obs_taxon, levels)
            if debug: print tip.Name, obs_taxon, exp_taxon
            recoveries.append(get_lineage_recovery(obs_taxon, exp_taxon))
        return recoveries, tipnames



    def _calc_tip_taxon_without_self_fast(self, tip):
        """return the taxon predicted for the tip.

        #zero tip.TaxArray, return,  then restore it.
        """
        ori_arrays = [tip.TaxArray.copy()]

        #nullify tip and recalc ancs
        tip.TaxArray[:] = 0
        ancestors = list(tip.ancestors())
        for anc in ancestors:
            ori_arrays.append(anc.TaxArray.copy())
            anc.TaxArray = self._get_tax_array_fitch(anc)
        if self._fitch_back:
            for anc in reversed(ancestors[:-1]): #skip root
                self._adjust_states(anc)
        result = self._calc_tip_taxon(tip)

        #restore .TaxArray of self and ancs
        tip.TaxArray = ori_arrays.pop()
        for anc in ancestors:
            anc.TaxArray = ori_arrays.pop()

        return result



















    ### obseleted below
    
    def _calcScoresByLeftOneOut(self, exclude_leftout_tip=False,):

        """return scores, tipnames for each tip with taxon.
        
        exclude_leftout_tip: if True, nullify the leftout tip.TaxArray
        """
        #predict all the taxa
        if not exclude_leftout_tip:
            taxa_predict = self.calcAllTipTaxa(exclude_taxa_tips=False)
        else:
            taxa_predict = self.calcAllTipTaxa(exclude_taxa_tips=False,
                    tiptaxon_getter = self._calc_tip_taxon_without_self)

        scores, tipnames = [], []
        for tip in self._root.iterTips():
            tipname = tip.Name
            tiptaxon = self._taxas.get(tipname, [])
            if not tiptaxon: continue #only score for known taxa

            #pretreat the tiptaxon for fair comparison
            tiptaxon = map(str.lower, tiptaxon)[:maxitems]
            #ignore subtax names
            tiptaxon = [tax for tax in tiptaxon
                    if not tax in self._rdp.SUB_TAXS]
            #truncate at tax names meaning uncertain
            for i, tax in enumerate(tiptaxon):
                if tax.startswith('incertae sedis')\
                        or tax.startswith('unclassified_'):
                    tiptaxon = tiptaxon[:i]
                    break

            #copare and score
            taxon_predict = taxa_predict[tipname]
            scores.append(self._score_taxon(taxon_predict, tiptaxon))
            #for debuging (tipname, obs, exp)
            tipnames.append((tipname, taxon_predict, tiptaxon))

        return scores, tipnames


    @staticmethod
    def _score_taxon(obs, exp, maxitems=None):
        """return a score (float w\in [-1, 1]), after comparing obs and exp.

        a ratio if right, a negative ratio if wrong.
        Note: no penalty for a predict longer than exp.
        Note: obs, exp should support a[:i]
        """
        if maxitems:
            exp = exp[:int(maxitems)]
        if not exp: return nan
        if not obs: return 0

        len_exp = len(exp)
        pairs = map(None, obs[:len_exp], exp)
        for same_count, (obs_item, exp_item) in enumerate(pairs):
            if obs_item != exp_item: #obs item is wrong or None
                break
        else: #all the same within len(exp)
            return 1

        #figure out the score
        result = same_count/len_exp
        if obs_item: #obs item unmatch the exp one.
            return result - 1
        else: #obs with less items
            return result
























##################
# deprecated below
from operator import itemgetter
from itertools import count
from py_util.misc import lace #zip after filling in the shorter seq
import re
multiwhite = re.compile('\s{3,}') #used to split the arb nds lines


def _first_two_words(line):
    """return the first two words or two Nones."""
    result = line.split()[:2]
    if len(result) != 2:
        return None, None
    return result
## moved to util.py
class RdpTaxonomy_old(object):
    """Provide methods for a RDP taxonomy.

    Implementation Notes:
    each {tax_array} is actually a multi bitsets with borders for each level (phylum, class, etc)
    
    Todo: maxlevel to treenode.breath_first?
    todo: maxlevel default to be nolimit.
    Todo: use sparse array or frozen sets for TaxArray.
    """
    #taxonomy levels
    Levels = dict(zip(
        ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
        count(1)))
    #tax names to be ignored: suborder subclass subsection
    SUB_TAXS = dict.fromkeys([
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
        'glycomycineae'])

    def __init__(self, tree, maxlevel=None, index_tree=True):
        """load the tree as .Root"""
        if not isinstance(tree, TreeNode):
            tree = DndParser(tree, TreeNode)
        self._root = tree

        if maxlevel is None:
            maxlevel = self.Levels['family']

        if index_tree:
            # _nodes_bf: node list in breath_first order
            self._nodes_bf, self._level_slices = \
                    self._index_tree(self._root, maxlevel)

            num_nodes = len(self._nodes_bf)
            self._empty_tax_array = zeros(num_nodes, int)


    @classmethod
    def fromTaxLines(cls, taxa_tree_lines, maxlevel=None, parser=_first_two_words):
        """construct taxtree from the text format.
        
        todo: fix maxlevel, no easy.
        """
        def _adjust_length_of_active_nodes(level):
            """fill or trim the active node list, if level != list.length"""
            delta = level - len(active_nodes)
            if delta > 0:
                #fill-in with a lineage of empty nodes
                for i in range(delta):
                    new_node = TreeNode('')
                    active_nodes[-1].addChild(new_node)
                    active_nodes.append(new_node)
            elif delta < 0:
                for i in range(-delta):
                    active_nodes.pop()
            #else: pass when delta == 0

        levels_get = cls.Levels.get #speed hack
        root = TreeNode()
        active_nodes = [root] #keep active nodes in level order
        for line in taxa_tree_lines:
            title, name = parser(line)
            level = levels_get(title)
            if not level:# or level > maxlevel: #not valid level title
                continue
            #if level > maxlevel:
            #    _adjust_length_of_active_nodes(maxlevel)
            #    continue

            _adjust_length_of_active_nodes(level)
            #add a new child
            new_node = TreeNode(name)
            active_nodes[-1].addChild(new_node)
            active_nodes.append(new_node)

        return cls(root, maxlevel)

    def _index_tree(self, root, maxlevel):
        """index the tree in breadthfirst order.

        return a list of nodes in breadfirst order, a list of slices for each
        level.
        todo: what if maxlevel if unlimited.
        todo: why need to set node._index node._level and modify node.Name?
        todo: why not just level_borders
        todo: will {for each level in bf_levels: for each node in level: do something}  easier.
        """
        ## node._index: index of the node in nodes_bf
        nodes_bf = []
        level_slices = []
        prev_i, prev_level = 0, 0
        for i, (node, level) in enumerate(self._root.breadth_first_with_level()):
            if level > maxlevel: break
            if level > prev_level:
                level_slices.append(slice(prev_i, i))
                prev_i, prev_level = i, level
            node._index, node._level = i, level
            #lowercase for non case sensitive search
            if node.Name:
                node.Name = node.Name.lower()
            nodes_bf.append(node)

        level_slices.append(slice(prev_i, None))
        return nodes_bf, level_slices

    def taxArrayFromTaxon(self, taxon):
        """return an array with 0/1 representing existence of a node.

        - taxon: a seq of tax names.

        Note: ignore subclasses and suborders and subsections
        """
        result = self._empty_tax_array.copy() #init as zeros
        if not taxon: #taxon as []
            return result
        taxon = map(str.lower, taxon)

        curr_nodes = self._root.Children #first level
        for tax in taxon:
            if tax in self.SUB_TAXS: continue #ignore the sub_ranks
            node = self._find_tax(tax, curr_nodes)
            if not node: #node not found with the name tax
                break
            #node found
            result[node._index] = 1
            curr_nodes = node.Children #next level, if found
        return result

    @staticmethod
    def _find_tax(tax, curr_nodes):
        """return the first node with .Name eq tax, or None

        Warn: if there are duplicated names in curr_nodes, only return the fst.
        """
        #tax = tax.lower() #should be done in outer functions for efficiency
        for node in curr_nodes:
            if node.Name == tax:
                return node
        #return None

    def taxonFromTaxArray(self, tax_array):
        """return a taxon from a tax array, or []."""
        nodes_bf = self._nodes_bf
        level_slices = self._level_slices[1:] #donot want the root
        tax_array = asarray(tax_array, int) #for unit test

        result = []
        for level_slice in level_slices:
            level_taxs = tax_array[level_slice]
            nz_idxs = level_taxs.nonzero()[0]
            if len(nz_idxs) != 1: #no tax or multiple taxs
                break
            node_idx = level_slice.start + nz_idxs[0]
            result.append(nodes_bf[node_idx].Name)
        return result

        
class TaxaTree(object):
    """Provide methods to init and prediction taxa based on a tree.
    """
    def __init__(self, tree=None, taxas=None):
        """init tree.
        """
        if tree:
            self.loadTree(tree)
        if taxas:
            self.loadTaxas(taxas)

    def loadTree(self, tree):
        """load the tree as .Root and synchronize tree with taxas if possible.
        """
        if not isinstance(tree, PhyloNode):
            tree = DndParser(tree)
        self._root = tree

        if hasattr(self, '_taxas'):
            self._initTreeTaxas()

    def loadTaxas(self, taxas):
        """load the taxas as .Taxas and synchronize tree with taxas if possible.
        """
        if not isinstance(taxas, dict):
            taxas = self._get_taxas_from_lines(taxas)
        self._taxas = taxas

        if hasattr(self, '_root'):
            self._initTreeTaxas()

    @staticmethod
    def _get_taxas_from_lines(name_taxa_lines):
        """return a dict of {tipname: taxa as a list}"""
        table = Table(name_taxa_lines, header=False)
        taxas = dict((name, taxa.split('; '))
                for name, taxa in table
                if taxa) #not include in result if no valid taxa
        return taxas

    def _initTreeTaxas(self):
        """synchronize tree and taxas.
        """
        taxas = self._taxas

        num_taxas = len(taxas)
        tipnames = taxas.keys()
        self._taxa_vectors = array(taxas.values(), object)

        taxa_defarrays = dict(zip(tipnames, identity(num_taxas, bool)))
        self._taxa_blankarray = taxa_blankarray = zeros(num_taxas, bool)

        #initiate .TaxaArray for tips from .Name
        for tip in self._root.iterTips():
            tip.TaxaArray = taxa_defarrays.get(tip.Name, taxa_blankarray)
            tip.Taxa = taxas.get(tip.Name, None)

        self._setNontipTaxaArrays()

    def _setNontipTaxaArrays(self):
        """refresh .TaxaArray for nontips from that of tips"""
        for node in self._root.iterNontips(mode='postorder'):
            node.TaxaArray = or_reduce([c.TaxaArray for c in node.Children])

    def calcTipTaxa(self, tip, **kw):
        """return the predicted taxa for a tip, using .TaxaArray of nontips.
        """
        for anc in tip.iterAncestors():
            curr_taxa = self._calcNontipTaxa(anc, **kw)
            if curr_taxa:
                return curr_taxa

    def _calcNontipTaxa(self, nontip, mode='common', **kw):
        #min_taxas=3, threshold=0.9, min_taxas_most=10):
        """
        mode='common': ['common', 'most', 'both']
        """
        valid_kw = ['min_taxas', 'min_taxas_most', 'threshold']
        assert set(kw.keys() + valid_kw) == set(valid_kw)
        if mode == 'common':
            return self._calcNontipTaxaCommon(nontip, **kw)
        elif mode == 'most':
            return self._calcNontipTaxaMost(nontip, **kw)
        elif mode == 'both':
            result = self._calcNontipTaxaCommon(nontip, **kw)
            if result:
                return result
            return self._calcNontipTaxaMost(nontip, **kw)
        else:
            raise ValueError()

    def _calcNontipTaxaCommon(self, nontip, min_taxas=3, **kw):
        """return the predicted taxa for a nontip node.
        
        min_taxas: min num of taxas to predict
        **kw: place holder, never used.
        """
        if sum(nontip.TaxaArray) < min_taxas:
            return None

        result = []
        taxas_to_use = self._taxa_vectors[nontip.TaxaArray]
        for samelevel_taxons in zip(*taxas_to_use):
            if len(set(samelevel_taxons)) > 1:
                break
            result.append(samelevel_taxons[0])
        return result

    def _calcNontipTaxaMost(self, nontip,
            threshold=0.9, min_taxas_most=10, **kw):
        """return the predicted taxa for a nontip node.
        
        min_taxas_most=10: min num of taxas to predict
        threshold=1: min majority needed to predict.
        **kw: place holder, never used.
        """
        num_taxas = sum(nontip.TaxaArray)
        if num_taxas < min_taxas_most:
            return None

        result = []
        scores = []
        taxas_to_use = self._taxa_vectors[nontip.TaxaArray]
        for samelevel_taxs in lace(*taxas_to_use):
            freqs = count_freqs(samelevel_taxs)
            #if len(freqs) == 1:
            #    result.append(freqs.keys()[0])
            #    scores.append(1)
            #elif threshold != 1:
            max_tax, max_freq = sorted(freqs.items(), key=itemgetter(-1))[-1]
            if max_tax and max_freq/num_taxas >= threshold:
                result.append(max_tax)
                scores.append(max_freq/num_taxas)
            else:
                break
        return result

    #for efficiency O(n)
    def calcAllTipTaxas(self, exclude_taxa_tips=True, **kw):
        """return a dict of {notaxa-tipname: taxa-predict}"""
        #self._setNontipTaxaArrays()
        self._setNontipTaxas(**kw)
        result = {}
        for tip in self._root.iterTips():
            if exclude_taxa_tips and tip.Taxa: continue
            result[tip.Name] = self._fetchTipTaxa(tip)
        return result

    def _setNontipTaxas(self, **kw):
        """set the .Taxa for all nontip nodes, using .TaxaArray."""
        for node in self._root.iterNontips():
            node.Taxa = self._calcNontipTaxa(node, **kw)

    def _fetchTipTaxa(self, tip):
        """return the first valid .Taxa from the tip's ancestors, or None.
        """
        for anc in tip.iterAncestors():
            if anc.Taxa:
                return anc.Taxa

    def _calcScoresByLeftOneOut(self, exclude_leftout_tip=False, **kw):
        """return scores, tipnames for each tip."""
        scores = []
        tip_taxas = [(t.Name, t.Taxa) for t in self._root.iterTips()]

        if not exclude_leftout_tip:
            taxas_predict = self.calcAllTipTaxas(exclude_taxa_tips=False, **kw)
        else:
            taxas_predict = {}
            for tip in self._root.iterTips():
                #remove tip taxa from .TaxaArray of ancestors
                for anc in tip.iterAncestors():
                    anc.TaxaArray -= tip.TaxaArray
                taxas_predict[tip.Name] = self.calcTipTaxa(tip, **kw)
                #restore .TaxaArray of ancestors
                for anc in tip.iterAncestors():
                    anc.TaxaArray += tip.TaxaArray

        scores, tipnames = [], []
        for tipname, tiptaxa in tip_taxas:
            if not tiptaxa: continue
            tipnames.append(tipname)
            try: taxa_predict = taxas_predict[tipname]
            except: set_trace()
            if taxa_predict:
                scores.append(_score_taxon(taxa_predict, tiptaxa))
            else:
                scores.append(0)
        return scores, tipnames

def _score_taxon(obs, exp, maxitems=None, filter_=None, map_=str.lower):
    """return a score (float), after comparing obs and exp.

    a ratio if right, a negative ratio if wrong or too long.
    Note: obs, exp expected to be iterables.

    todo: remove maxitems, filter_ and map_??
    """
    if maxitems:
        maxitems = int(maxitems)
        obs, exp = obs[:maxitems], exp[:maxitems]
    if filter_:
        obs = filter(filter_, obs)
        exp = filter(filter_, exp)
    if map_:
        obs = map(map_, obs)
        exp = map(map_, exp)

    pairs = lace(obs, exp)
    num_pairs = len(pairs)

    for same_count, pair in enumerate(pairs):
        if pair[0] != pair[1]:
            break
    else: #all the same
        return 1

    raw_score = same_count/num_pairs
    #return raw_score if (obs[same_count] is None) else raw_score-1
    if pair[0] is None: #obs shorter
        return raw_score
    elif pair[1] is None: #obs is incorrectly longer
        return raw_score - 1
    else: #obs item unmatch the exp one.
        return raw_score - 1


























#################################
# deprecated here below
class TaxaNode(PhyloNode):
    """PhyloNode with two extra attrs: .Taxa, .Acc."""
    def __init__(self, *a, **kw):
        super(TaxaNode, self).__init__(*a, **kw)
        self.Taxa = None
        self.TaxaArray = None
    #tobe removed here below for this class
        self.Acc = None
        self._tip_taxas = None #cache for efficency, tobe moved to Tree

    @staticmethod
    def _has_taxa(node): return node.Taxa

    @staticmethod
    def _no_taxa(node): return not node.Taxa

    def iterTaxaTips(self):
        """return each tip with .Taxa info."""
        return self.iterTips(filter_func=self._has_taxa)

    def iterNontaxaTips(self):
        """return each tip without .Taxa info."""
        return self.iterTips(filter_func=self._no_taxa)

    def tipTaxas(self, refresh=False):
        """return a list of taxas from tips."""
        if refresh or not self._tip_taxas:
            self._tip_taxas = [n.Taxa for n in self.iterTaxaTips()]
        return self._tip_taxas


def tree_setKnownTaxas(tree, tipname_taxas):
    """set the .Taxa for tips using the lookup dict.

    tree: a TaxaNode instance.
    name_taxas: a dict of {tipname: taxa as list}
    """
    for tip in tree.iterTips():
        tip.Taxa = tipname_taxas.get(tip.Name)

def tree_getUnknownTaxas(tree, min_taxas=9, threshold=None, return_more=False):
    """return a dict of {tipname: (num_taxas, [common_taxs], [(tax,
    percent),])}

    tree: a TaxaNode instance.
    threshold=None: a ratio threshold which will be count as tax.
    """
    result = {}
    for tip in tree.iterNontaxaTips():
        #node, parent = tip, tip.Parent
        #while parent is not None:
        for anc in tip.iterAncestors():
            taxas = anc.tipTaxas()
            if len(taxas) >= min_taxas:
                break
        else: #no .Taxa available
            result[tip.Name] = (0, [], [])
            continue

        num_taxas = len(taxas)
        common_taxs, max_taxs = [], []
        for samelevel_taxs in lace(*taxas):
            freqs = count_freqs(samelevel_taxs)
            if len(freqs) == 1:
                common_taxs.append(freqs.keys()[0])
            else:
                tax, max_freq = sorted(freqs.items(), key=itemgetter(-1))[-1]
                if threshold and max_freq/num_taxas >= threshold:
                    max_taxs.append((tax, max_freq))
        result[tip.Name] = (num_taxas, common_taxs, max_taxs)

    if not return_more:
        #only return the taxs
        result = dict((k, v[1]+ [t for t, f in v[2]])
                for k, v in result.items())

    return result

        
# but it is not mapped to RDP system yet??
def get_tipname_taxas(nds_lines, tab=True):
    """return a dict of {tipname: taxa as list} from Arb nds lines.
    """
    if tab:
        table = Table(nds_lines, header=False)
        result = dict((name, taxa.split('; '))
                for name, taxa in table)
        return result

    #if space delimited
    result = {}
    for line in nds_lines:
        #each line is in the format " name     tax1; tax2; tax3   "
        name, taxa_str, empty = multiwhite.split(line)
        name = name.strip() #there is a space before it
        taxa = taxa_str.split('; ')
        result[name] = taxa
    return result

def eval_taxa_result(tree, tipname_taxas, del_ratio=0.2, 
        func=tree_getUnknownTaxas, **kw):
    """return a score matrx of [[tipname, score, obs, exp],]

    del_ratio: the ration of items to delete from tipname_taxas
    """
    items = array(tipname_taxas.items(), object)
    num_items = len(items)
    num_dels = int(num_items * del_ratio)

    #rand remove certain number of items to test
    keep_idxs = permutation(num_items)[:-num_dels]
    test_mapping = dict(items[keep_idxs])

    tree_setKnownTaxas(tree, test_mapping)
    obs_mapping = func(tree, **kw)

    result = []
    for k, obs_taxa in obs_mapping.items():
        try:
            ori_taxa = tipname_taxas[k]
        except KeyError: #when k has no ori_taxa
            continue
        score = _score_taxon(obs_taxa, ori_taxa)
        result.append([k, score, obs_taxa, ori_taxa])

    return result




