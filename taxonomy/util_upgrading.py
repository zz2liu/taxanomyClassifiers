"""util.py

Similar Codes:
    lineage_tree.py (newer)

02/02/08: mv calc_taxa_recovery here from taxonomy_from_tipnames.py

todo: it is better to keep a normal tree from the taxonomy, with rank
infomation.  and return the lineage with ranks.
-> remove missing_word option from tree_search later.

todo: move .fromTaxLines to parse.tree
todo: Decide the behavior on the following cases
    1. subclass and suborder; -> make sub-tax-names for bergy7.8
    2. missing order and family;
    3. missing after class
    AY193128|S000401279 Root;Bacteria;OP11;OP11_genera_incertae_sedis
todo make ._index_tree a general tree function.
todo: make a Taxonomy(TreeNode)

"""
from pdb import set_trace
from itertools import count
from numpy import asarray, zeros
from py_util.tree_new import (DndParser, TreeNode, search, lineage,
        breadth_first, traverse)
from py_util.dict_ import calc_freqs

TITLES = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
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

def index_tree_by_level(root):
    """index the tree in breadthfirst order, assign ._index, ._level to each
    node.  return (nodes in breadfirst order, slices by level).

    todo: what if maxlevel if unlimited.
    todo: why need to set node._index node._level and modify node.Name?
    todo: why not just level_borders
    todo: will {for each level in bf_levels: for each node in level: do something}  easier.
    """
    ## node._index: index of the node in nodes_bf
    nodes_bf = []
    level_slices = []
    prev_i, prev_level = 0, 0
    for i, (node, level) in enumerate(breadth_first(root, with_level=True)):
        #if level > maxlevel: break  #quickfix while run_fitch.py
        if level > prev_level:
            level_slices.append(slice(prev_i, i))
            prev_i, prev_level = i, level
        node._index, node._level = i, level

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

def _first_two_words(line):
    """return the first two words or two Nones."""
    result = line.split()[:2]
    if len(result) != 2:
        return None, None
    return result

#todo make all methods to be class methods.
class Taxonomy(TreeNode):
    Levels = dict(zip(TITLES, count(1)))
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

    Implementation Notes:
    each {tax_array} is actually a multi bitsets with borders for each level (phylum, class, etc)

    Fixit: maxlevel doesn't work for bergy7.8??
    Todo: maxlevel to treenode.breath_first?
    todo: maxlevel default to be nolimit.
    Todo: use sparse array or frozen sets for TaxArray.
    """
    #taxonomy levels
    Levels = dict(zip(TITLES, count(1)))

    #tax names to be ignored: suborder subclass subsection
    SUB_TAXS = SUB_TAX_NAMES

    def __init__(self, tree, maxlevel='family', index_tree=True,
            verbose=False):
        """load the tree as .Root

        - maxlevel='family'
        """
        self.Verbose = verbose
        if not isinstance(tree, TreeNode):
            tree = DndParser(tree, TreeNode)
        self._root = tree

        if isinstance(maxlevel, str):
            maxlevel = self.Levels[maxlevel.lower()]

        if index_tree:
            # _nodes_bf: node list in breath_first order
            self._nodes_bf, self._level_slices = \
                    self._index_tree(self._root, maxlevel)

            num_nodes = len(self._nodes_bf)
            self._empty_tax_array = zeros(num_nodes, int)

    # used in parse_rdp.py
    @classmethod
    def fromRawTaxTree(cls, raw_tree, **kw):
        """construct taxtree from a raw tax tree with .params['rank'].

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
        return cls(raw_tree, **kw)
            
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

    def _index_tree_old(self, root, maxlevel):
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
            #if level > maxlevel: break  #quickfix while run_fitch.py
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

    def _index_tree(self, lower=True, maxlevel='genus'): #maxlevel never used
        """index the tree in breadthfirst order.
        """
        root = self._root
        if lower:
            for n in traverse(root):
                if not n.Name: continue
                n.Name = n.Name.lower()
        nodes, level_slices = index_tree_by_level(root)
        return nodes, level_slices
        
    def taxArrayFromTaxon(self, taxon):
        """return an array with 0/1 representing existence of a node.

        - taxon: a seq of tax names.

        assuming ._index have been assigned to the taxonomy nodes.
        Note: ignore subclasses and suborders and subsections
        """
        result = self._empty_tax_array.copy() #init as zeros
        if not taxon: #taxon as []
            return result
        taxon = [tax and tax.lower() or '' for tax in taxon]

        curr_nodes = self._root.Children #first level
        for tax in taxon:
            if tax in self.SUB_TAXS: continue #ignore the sub_ranks
            node = self._find_tax(tax, curr_nodes)
            if node is None: #fix a bug, node not found with the name tax
                break
            #node found
            try:result[node._index] = 1
            except: set_trace()
            curr_nodes = node.Children #next level, if found
        return result


    @staticmethod
    def _find_tax(tax, curr_nodes):
        """return the first node with .Name eq tax, or None

        Warn: if there are duplicated names in curr_nodes, only return the fst.
        """
        for node in curr_nodes:
            if node.Name == tax:
                return node


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

        

def calc_taxa_recovery(obs_taxa, exp_taxa, num_levels=5, 
        strict=False, verbose=False):
    """return a list of recovery rates for each level.

    obs_taxa, exp_taxa: expected to be dict of {seqname: taxon as a list of strs}
    """
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


