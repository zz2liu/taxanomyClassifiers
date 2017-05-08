"""util.py

Similar Codes:
    lineage_tree.build_genus_tree

12/02/07: add a Taxonomy(TreeNode)
02/02/08: mv calc_taxa_recovery here from taxonomy_from_tipnames.py
02/14/08 deprecate maxlevel in RdpTaxonomy; deprecate calc_taxa_recovery
02/14/08 seperate several methods to be functions.
02/14/08 do not lower the node names while _index_tree now.
02/14/08 : mv SUBTAXES to data
02/15/08 : tree_from_taxlines generalized to parse.tree

todo: merge Taxonomy and RdpTaxonomy?  How to store the _level_slices and
    bf_nodes, if RdpTaxomy becomes a TreeNode?
todo: make taxarray as bool array or bit array (later) or sparse array or
    frozenset?  The reason i kept it as int array if for future count nodes.
todo: return level starts in ._index_tree and use array_split in taxonFromArray
todo: unify the node searching of TreeLCA and TreeFitch
todo: it is better to keep a normal tree from the taxonomy, with rank
infomation.  and return the lineage with ranks.
-> remove missing_word option from tree_search later.

todo: Decide the behavior on the following cases
    1. subclass and suborder; -> make sub-tax-names for bergy7.8
    2. missing order and family;
    3. missing after class
    AY193128|S000401279 Root;Bacteria;OP11;OP11_genera_incertae_sedis

"""
from functools import partial
from pdb import set_trace
from itertools import count
from warnings import warn
from numpy import asarray, zeros

from py_util.tree_new import (DndParser, TreeNode, search, search_lite,
        lineage, traverse)
from py_util.parse.tree import tree_from_name_level_lines
from py_util.dict_ import calc_freqs
from data import SUB_TAX_NAMES

breadth_first_with_level = partial(traverse,
        mode='breadth_first', with_level=True)
TITLES = ['domain', 'phylum', 'class', 'order', 'family', 'genus']#, 'species']
LEVELS = dict(zip(TITLES, count(1)))


###
# parse tree from RDP indented tax lines
_levels_get = LEVELS.get #-> 1..7 or None
def _first_two_words(line):
    """return the first two words or two Nones."""
    result = line.split()[:2]
    if len(result) != 2:
        return None, None
    return result

def _parse_rdp_tax_line(line):
    title, name = _first_two_words(line)
    level = _levels_get(title)
    return name, level

tree_from_taxlines = partial(tree_from_name_level_lines,
        parser=_parse_rdp_tax_line)


def standize_ranked_tree(raw_tree):
    """modify the raw tree by padding or removing nodes.
    Assume .params['rank'], .iterNodes, .removeSelf for each node.

    Make the tree strictly with levels DPCOFG...
    Note: the raw tree should have a method .removeSelf
    """
    #remove the subclass and suborder nodes (with a rank 3.5 or 4.5
    nodes_to_remove = [n for n in raw_tree.iterNodes(include_self=False)
            if n.params['rank'] % 1]
    for node in nodes_to_remove:
        node.removeSelf()
    #add nodes to gaps
    nodes_to_padabove = [n for n in raw_tree.iterNodes(include_self=False)
            if n.params['rank'] - n.Parent.params['rank'] > 1]
    for node in nodes_to_padabove:
        delta = n.params['rank'] - n.Parent.params['rank']
        assert delta > 1 and delta%1 == 0
        parent = node.Parent
        #fill-in with a lineage of empty nodes
        for i in range(delta-1):
            new_node = TreeNode('')
            new_node.Parent = parent
            new_node.params['rank'] = -1
            parent = new_node
        node.Parent = parent
    return raw_tree

def index_tree_by_level(root):
    """index tree in breadth_first order, setting ._index
    Return (nodes, slices of levels) include root.
    """
    ## node._index: index of the node in nodes_bf
    nodes_bf = []
    level_slices = []
    prev_i, prev_level = 0, 0
    for i, (node, level) in enumerate(breadth_first_with_level(root)):
        #if level > maxlevel: break  #quickfix while run_fitch.py
        if level > prev_level:
            level_slices.append(slice(prev_i, i))
            prev_i, prev_level = i, level
        node._index = i
        nodes_bf.append(node)
    level_slices.append(slice(prev_i, None))
    return nodes_bf, level_slices

def get_single_genus_ids(dpcofg, i=-1):
    """return the keys with a single value[i].

    - dpcofg: a {key, value as a list with same lengths}
    - i: the index of value to count freqs.
    """
    genera = [v[i] for v in dpcofg.values()]
    genus_freqs = calc_freqs(genera)
    result = [id for id, lin in dpcofg.items()
            if genus_freqs[lin[i]] == 1]
    return result

#todo make all methods to be class methods.
class Taxonomy(TreeNode):
    Levels = LEVELS
    def taxNodeFromTaxon(self, taxon, taxon_with_root=True,
            ignore_case=True, include_subranks=False,
            verbose=False):
        """return a node with a list of tax names.

        Note: ignore subclasses and suborders and subsections
        """
        if ignore_case: #must before searching subranks
            taxon = map(str.lower, taxon)
        taxon = [tax or '' for tax in taxon]
        if not include_subranks:
            taxon = [tax for tax in taxon if tax not in SUB_TAX_NAMES]

        result = search(self, taxon, ignore_case=ignore_case,
                verbose=verbose, words_with_root=taxon_with_root,
                missing_word='') #to find the genus of Cyanobacteria 
                        #misssing order, family
        return result

    def getTaxon(self, taxon_with_root=True):
        """return a taxon from a tax array, or []."""
        result = [a.Name for a in lineage(self)]
        if not taxon_with_root:
            return result[1:]
        else:
            return result


class RdpTaxonomy(object):
    """Provide methods for a RDP taxonomy.
    """
    #taxonomy levels
    Levels = LEVELS

    #tax names to be ignored: suborder subclass subsection
    SUB_TAXS = SUB_TAX_NAMES

    def __init__(self, tree, **kw):
        """load the tree as .Root

        **kw: (deprecated)
            - maxlevel='family'
            - index_tree: only useful for test.
        """
        if 'maxlevel' in kw:
            warn('do not use maxlevel', DeprecationWarning)
        maxlevel = kw.pop('maxlevel', 'family')
        if isinstance(maxlevel, str):
            maxlevel = self.Levels[maxlevel.lower()]
        index_tree = kw.pop('index_tree', True)
        assert not kw

        if not isinstance(tree, TreeNode):
            tree = DndParser(tree, TreeNode)
        self._root = tree

        if index_tree:
            # _nodes_bf: node list in breath_first order
            self._nodes_bf, self._level_slices = \
                    self._index_tree()

            num_nodes = len(self._nodes_bf)
            self._empty_tax_array = zeros(num_nodes, int)

    def _index_tree(self, **kw):
        """index the tree in breadthfirst order.
        return (nodes in breadfirst order, slices for levels)
        """
        if 'maxlevel' in kw:
            warn('', DeprecationWarning)
        kw.pop('maxlevel', None)
        assert not kw

        root = self._root
        return index_tree_by_level(root)

    def taxArrayFromTaxon(self, taxon):
        """return an array with 0/1 representing existence of a node.

        - taxon: a seq of tax names.

        assuming ._index have been assigned to the taxonomy nodes.
        Note: ignore subclasses and suborders and subsections
        """
        result = self._empty_tax_array.copy() #init as zeros
        if not taxon: #taxon as []
            return result
        taxon = [tax and tax.lower() or '' #None -> ''
                for tax in taxon
                if not tax in self.SUB_TAXS]

        nodes = search_lite(self._root, taxon, lower=True)
        for n in nodes:
            result[n._index] = 1
        return result

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





    ###
    # deprecatec methods
    # used in parse_rdp.py
    @classmethod
    def fromRawTaxTree(cls, raw_tree, **kw):
        """construct taxtree from a raw tax tree with .params['rank'].
        Make the tree strictly with levels DPCOFG...
        """
        standize_ranked_tree(raw_tree)
        return cls(raw_tree, **kw)
            
    @classmethod
    def fromTaxLines(cls, taxa_tree_lines, maxlevel=None, **kw):
        """construct taxtree from the text format.
        
        todo: fix maxlevel, no easy.
        """
        if maxlevel is not None:
            warn('maxlevel', DeprecationWarning)

        root = tree_from_taxlines(taxa_tree_lines, **kw)
        return cls(root, maxlevel=maxlevel)














####
# deprecated
def calc_taxa_recovery(obs_taxa, exp_taxa, num_levels=5, 
        strict=False, verbose=False):
    """return a list of recovery rates for each level.

    obs_taxa, exp_taxa: expected to be dict of {seqname: taxon as a list of strs}
    """
    warn('use dpcofg_table_to_recovery.py', DeprecationWarning)
    if isinstance(obs_taxa, dict):
        #no penalty for those taxon not observed (clusted to 'Unclassified')
        common_keys = set(obs_taxa.keys()) & set(exp_taxa.keys())
        taxon_pairs = [(obs_taxa[k], exp_taxa[k]) for k in common_keys]
    else: #taxa as vector
        taxon_pairs = zip(obs_taxa, exp_taxa)

    #init the score_matrix with NotKnowns (nan)
    num_pairs = len(taxon_pairs)
    score_matrix = zeros((num_pairs, num_levels), float)
    score_matrix.fill(nan)

    #fill the score_matrix:
    #where e is known: if o == e, 1; if o!= e, 0
    for i, (obs, exp) in enumerate(taxon_pairs):
        #no penalty for longer obs list.
        exp = exp[:num_levels]
        obs = obs[:len(exp)]
        obs += [None] * (len(exp)-len(obs))

        after_error = False
        for j, (o, e) in enumerate(zip(obs, exp)):
            if after_error:
                score_matrix[i, j] = 0
                continue
            if o == e:
                score_matrix[i, j] = 1
            else: #o != e or o is None
                score_matrix[i, j] = 0
                after_error = True
    equals = (score_matrix == 1).sum(0) #for each level
    unknowns = isnan(score_matrix).sum(0)
    percentages = equals/(num_pairs - unknowns)
    if verbose:
        print 'equals: %s, unknowns: %s' % (equals, unknowns)
    return percentages


