"""taxonomy_tree.py

Similar Code??
    taxonomy_from_rdp_classifier.py

"""
from numpy import asarray
from py_util.tree_new import DndParser

def _first_two_words(line):
    """return the first two words or two Nones."""
    result = line.split()[:2]
    if len(result) != 2:
        return None, None
    return result

class RdpTaxonomy(object):
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

        
