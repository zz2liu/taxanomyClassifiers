"""Provide TaxaTree to parse rdp taxonomy.

02/15/08 : looked into this code.

todo: replace TreeFromRdpLines with parse.tree.tree_from_name_level_lines, or
that from taxonomy_from_tipnames.

todo: the name mapping mv to data

todo: rewrite the _findTaxNode and taxFromTaxnode according to that in
taxonomy_from_tipnames.

Todo: instead of going up, map to a known RDP name when got a unknown groupname.
todo: move treeFromLines and tax recovery code to other modules.

#not standarized yet.
Terms used in this module:
    taxon: ['bacteria', 'firmicutes', 'mollicutes']
    taxa: a list of taxon.
    tax: one item in taxon. eg. 'mollicutes'
    taxs: a list of tax not forming a taxon. eg. ['mollicutes', 'bacilli']

"""
from pdb import set_trace
from itertools import chain, count as icount
import re
from warnings import warn
#from cogent.parse.tree import DndParser
from cogent.parse import tree_xml, newick
from pprint import pprint
from py_util.tree_new import DndParser, TreeNode
from py_util.dict_ import multival_dict

def _get_title_name(line):
    """return title and name from a RDP taxonomy line.
    
    - line: ex. '  class Mollicutes (0/1149/0)'
    """
    title_name = line.split(' (', 1)[0].strip()
    try:
        title, name = title_name.split(None, 1)
    except ValueError: #unclassified_Mycoplasmataceae (0/5/0)
        title, name = None, title_name
    return title, name

def TreeFromRdpLines(taxa_tree_lines, parser=_get_title_name):
    """return the root of a tree parsed from indented RDP taxonomy text.

    - taxa_tree_lines: file got from RDP website. head of it:
        domain Bacteria  (0/99183/0)  (selected/total/search matches)
                phylum Aquificae (0/558/0)
                    class Aquificae (0/558/0)
                        order Aquificales (0/558/0)
                            family Aquificaceae (0/529/0)
                                genus Aquifex (0/6/0)
    - parser: a function returning title and name from a line.

    Using: TITLES, TreeNode
    """
    #TITLES def in class TaxonomyFromGroupnames
    title_levels = dict(zip(TITLES, icount(1)))
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
    return root

def get_self_and_ancestor_names(node):
    """return a list of Node names (if any) of [self, ancestors].
    """
    return [n.Name
            for n in chain([node], node.ancestors())
            if n.Name]


_pace_to_rdp_map = {
        'Gamma proteobacteria A': 'Gammaproteobacteria',
        'Delta proteobacteria': 'Deltaproteobacteria',
        'Epsilon proteobacteria': 'Epsilonproteobacteria',
        'Hyphomicrobiaceae2': 'Hyphomicrobiaceae',
        'Gemmimonas': 'Gemmatimonas',
        'Haloanaerobiales':  'Halanaerobiales',
        'Haloanaerobiaceae':  'Halanaerobiaceae',
        'symbionts': 'Symbiotes',
        'Symbionts': 'Symbiotes',
        #with a -
        'Deinococcus-Thermus':  'Deinococcus-Thermus',

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
_subseps = re.compile('[ .-]') #move out of function for efficiency

def pacename_to_rdpname(name):
    """return the best RDP counterparts for a Pace tax name.

    Note: no warranty that the result exist in RDP.
    """
    #rename with a dict if possible
    try:
        return _pace_to_rdp_map[name]
    except KeyError:
        #rescue some names from the first word
        if ('/' not in name) and (',' not in name):
            name = _subseps.split(name)[0]
        return name



class TaxonomyFromGroupnames(object):
    """Provide taxonomy recall from the groupnames in phylogeny tree provided.
    """
    TITLES = 'domain phylum class order family genus species'.split()
    PACENAME_MODIFY_MAP = {#modify for consistency for duplicated names
        'Aquificae.2':  'Aquificae',
        'Aquificae':  'Aquificae-2', #only one class
        'Cyanobacteria.2': 'Cyanobacteria',
        'Cyanobacteria': 'Cyanobacteria-2', #only one class
        'Nitrospira.2': 'Nitrospira',
        'Nitrospira': 'Nitrospira-2', #only one class, genus will become class!
        'Acidiphilium.2': 'Acidiphilium',
        'Acidiphilium': 'Acidiphilium-2', #no duplicates
        'Sulfobacillus.2': 'Sulfobacillus',
        'Sulfobacillus': 'Sulfobacillus-2', #no duplicates
        'Synergistes.2': 'Synergistes', 
        'Synergistes': 'Synergistes-2', #no duplicates
        #save classes for monkey 454 data
        'Fibrobacteres': 'Fibrobacterales', #only class-order
        'Verrucomicrobia': 'Verrucomicrobiales', #only class-order
        'Chlorobi': 'Chlorobiales', #only class-order
        'Deferribacteres': 'Deferribacterales', #only class-order
        'Actinobacteria': 'Actinobacteria-', #only class
        }


    def __init__(self, tax_tree, phylo_tree):
        """initiate with the taxonomy tree and phylogeny tree.

        - tax_tree, phylo_tree expected to be subclass of TreeNode.
        """
        self._tax_tree = tax_tree
        #map each tax_node name to a list of tax_node (names not unique!)
        self._tax_nodes = multival_dict((n.Name, n)
                for n in self._tax_tree.breadth_first()
                if n.Name)

        #modify the node names, make nonnames empty.
        for node in phylo_tree.traverse():
            if not node.Name or node.Name.startswith('edge.'):
                node.Name = ''
        self._phylo_tree = phylo_tree

    def find_tax_node(self, name):
        """return a node found in the taxonomy tree by a pace name.
        
        - name: a group-name in Pace system.  Assuming that there is a '-' in
          lower level names if names duplicates.

        Dependants: pacename_to_rdpname
        """
        tax_nodes = self._tax_nodes

        name = name.strip("'")
        name = self.PACENAME_MODIFY_MAP.get(name, name)

        rdp_lookup = pacename_to_rdpname(name)
        nodes_found = tax_nodes.get(rdp_lookup, [None])
        if len(nodes_found) == 1:
            return nodes_found[0]
        else: #len(nodes_found) > 1:
            if '-' in name:
                return nodes_found[1]
            else:
                return nodes_found[0]
            if len(nodes_found) > 2:
                warn('the lowest levels of name %s is ignored' % repr(name))

    def _node_to_taxon(self, node, strict=False):
        """return group_names_from_root, group_name_mapped, taxon.
        """
        #try to map the closest group-name to a rdp node, if failed try
        #each upstream group-names
        group_names_to_root = get_self_and_ancestor_names(node)
        for name in group_names_to_root:
            tax_node = self.find_tax_node(name)
            if tax_node is not None: #TreeNode make overide bool!
                taxon = get_self_and_ancestor_names(tax_node)[::-1]
                return group_names_to_root[::-1], name, taxon
        else: #no group names from node up
            msg = 'No matched taxa_node for %s' % repr(node.Name)
            if strict:
                raise ValueError(msg)
            else:
                warn(msg)
                return group_names_to_root[::-1], '', ''

    def eachGroupnamesTaxon(self, verbose=False, strict=True):
        """yield (groupnames_from_root, mapped_name, taxon) for each groupname.
        """
        #for each nontip node with a group name
        for node in self._phylo_tree.traverse():
            if not node.Children or not node.Name:
                continue #skip leaves and nontax nodes
            group_names_from_root, name_mapped, taxon = \
                    self._node_to_taxon(node)
            yield group_names_from_root, name_mapped, taxon


    def eachTipTaxon(self, selector=None, strict=False, return_others=False):
        """yield (tipname, taxon) for each tip.

        - selector: a predicate(tip) -> bool, True to be selected.
        """
        for tip in self._phylo_tree.tips():
            if selector and not selector(tip): continue
            group_names_from_root, name_mapped, taxon = \
                    self._node_to_taxon(tip.Parent, strict=strict)
            if return_others:
                yield tip.Name, taxon, name_mapped, group_names_from_root
            else:
                yield tip.Name, taxon




TITLES = TaxonomyFromGroupnames.TITLES
