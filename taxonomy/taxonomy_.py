"""Provide TaxaTree to parse rdp taxonomy.
"""
from pdb import set_trace
from itertools import chain
import re
#from cogent.parse.tree import DndParser
from cogent.parse import tree_xml, newick
from pprint import pprint
from py_util.tree_new import DndParser, TreeNode

def _first_two_words(line):
    result = line.split()[:2]
    if len(result) != 2:
        return None, None
    return result

def _taxnames(node):
    """return a list of Node names from top level to node."""
    return [n.Name
            for n in reversed([node] + node.ancestors())
            if n.Name]#[1:]

_subseps = re.compile('[ .-]')

_pace_to_rdp = {
        'Gamma proteobacteria A': 'Gammaproteobacteria',
        'Hyphomicrobiaceae2': 'Hyphomicrobiaceae',
        #subclasses and suborders
        'Cystobacterineae': 'Myxococcales',
        'Sorangineae': 'Myxococcales',
        'Nannocystineae': 'Myxococcales',

        'Acidimicrobidae': 'Actinobacteria',
        'Acidimicrobineae': 'Acidimicrobiales',
        'Rubrobacteridae': 'Actinobacteria',
        'Rubrobacterineae': 'Rubrobacterales',
        'Coriobacteridae': 'Actinobacteria',
        'Coriobacterineae': 'Coriobacteriales',
        'Sphaerobacteridae': 'Actinobacteria',
        'Sphaerobacterineae': 'Sphaerobacterales',
        'Actinobacteridae': 'Actinobacteria',

        'Actinomycineae': 'Actinomycetales',
        'Micrococcineae': 'Actinomycetales',
        'Corynebacterineae': 'Actinomycetales',
        'Micromonosporineae': 'Actinomycetales',
        'Propionibacterineae': 'Actinomycetales',
        'Pseudonocardineae': 'Actinomycetales',
        'Streptomycineae': 'Actinomycetales',
        'Streptosporangineae': 'Actinomycetales',
        'Frankineae': 'Actinomycetales',
        'Glycomycineae': 'Actinomycetales',
        }

def _pace_rdp_mapping(name):
    #rename with a dict
    if name in _pace_to_rdp:
        return _pace_to_rdp[name]

    #rename those '..proteobacteria'
    if name.endswith('proteobacteria') and len(name.split()) == 2:
        return ''.join(name.split())

    #rescue some names from the first word
    if ('/' not in name) and (',' not in name):
        name = _subseps.split(name)[0]

    return name

def _get_taxa_node(name, taxa_nodes):
    taxa_node = taxa_nodes.get(name)

    #rescue with a mapping
    if taxa_node is None:
        name = _pace_rdp_mapping(name)
        taxa_node = taxa_nodes.get(name)
    return taxa_node

def _is_nonname(name):
    return (not name) or (name.startswith('edge.'))

class Taxonomy(object):
    Titles = 'domain phylum class order family genus species'.split()
    EmptyNode = TreeNode()
    def __init__(self, taxa_tree=None, phylo_tree=None):
        """taxa_tree, phylo_tree expected to be TreeNode."""
        self._taxa_tree = taxa_tree
        self._phylo_tree = phylo_tree
        #if taxa_tree:
        #    self.loadTaxaTree(taxa_tree, parser)
        #if phylo_tree:
        #    self.loadPhyloTree(phylo_tree)

    @property
    def TaxaTree(self): return self._taxa_tree

    @property
    def PhyloTree(self): return self._phylo_tree

    def loadTaxaTree(self, taxa_tree_lines, parser=_first_two_words):
        num_titles = len(self.Titles)
        title_levels = dict(zip(self.Titles, range(1, num_titles+1)))
        root = TreeNode()
        active_nodes = [root]

        def _update_active_nodes(level):
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
            #pass when delta == 0

        for line in taxa_tree_lines:
            title, name = parser(line)
            level = title_levels.get(title)
            if level: #skip those non-valid-titles
                _update_active_nodes(level)
                #add a new child
                new_node = TreeNode(name)
                active_nodes[-1].addChild(new_node)
                active_nodes.append(new_node)

        self._taxa_tree = root


    def loadPhyloTree(self, phylo_tree_lines, is_nonname=_is_nonname):
        """load and setup self._phylo_tree"""
        #parse the treelines
        result = DndParser(phylo_tree_lines)
        #modify the node names, make nonnames empty.
        for node in result.traverse():
            if is_nonname(node.Name): node.Name=''
        self._phylo_tree = result

    def calcGroupnameTaxa(self, verbose=False):
        """return a dict of {groupname: tax, ...}
        """
        #postorder here should be safer: when names duplicate, the higher order kept
        taxa_nodes = dict((n.Name, n)
                for n in self._taxa_tree.traverse(self_before=False, self_after=True)
                if n.Name)

        groupnames_taxa = {} #for convenience
        result = {}
        #for each tax node
        for node in self._phylo_tree.traverse():
            if not node.Children or not node.Name:
                continue #skip leaves and nontax nodes
            name0 = node.Name

            #prepare available names for mapping
            self_ancestor_names = [n.Name
                    for n in [node] + node.ancestors()
                    if n.Name]
            #try to map the first group name to a rdp node
            for name in self_ancestor_names:
                taxa_node = _get_taxa_node(name, taxa_nodes)
                if taxa_node is not None:
                    break
            else: raise ValueError('No matched taxa_node for %s' % repr(name0))

            result[name0] = _taxnames(taxa_node)
            groupnames_taxa[tuple(_taxnames(node))] = _taxnames(taxa_node)

        self._groupname_taxa = result
        return result, groupnames_taxa

    def calcSeqnameTaxa(self, phylo_tree=None, group_taxs=None):
        """return a dict of {seqname: tax, ...}

        using self._groupname_taxa if group_taxs not provided
        use self._phylo_tree if phylo_tree not provided.
        """
        if phylo_tree:
            self.loadPhyloTree(phylo_tree)
        if not group_taxs:
            group_taxs = self._groupname_taxa #speed hack

        _seqname_taxnames = {} #for convenience
        result = {}
        for tip in self._phylo_tree.tips():
            seqname = tip.Name
            tax_names = _taxnames(tip)[:-1] or ['']

            result[seqname] = group_taxs.get(tax_names[-1], [])
            _seqname_taxnames[seqname] = tax_names

        self._seqname_taxnames = _seqname_taxnames
        
        set_trace()
        return result
