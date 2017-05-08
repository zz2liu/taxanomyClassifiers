"""convenient wrappers for taxonomy classes."""
from cPickle import dump
from py_util.tree_new import DndParser
from data import RDP_TREE_PATH
from taxonomy.taxonomy_from_groupnames import (TaxonomyFromGroupnames)

def dump_tip_taxa_from_groupnames(phylo_tree_path, pickle_path,
        tax_tree_path=RDP_TREE_PATH, selector=None, verbose=False):
    """return a dict of {tip: t axon}

    Warning: hard coded.
    """
    taxa_tree = DndParser(file(tax_tree_path))
    phylo_tree = DndParser(file(phylo_tree_path))
        #'%s/454A_rick_w_bacteria.ntree' % ARB_DATA))
    taxonomy = TaxonomyFromGroupnames(taxa_tree, phylo_tree)

    tip_taxa = dict(taxonomy.eachTipTaxon(selector=selector))
    result = dict((k, map(str.lower, v))
            for k,v in tip_taxa.iteritems())
    dump(result, file(pickle_path, 'w'))
    return result
