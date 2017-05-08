"""taxononomy_from_nearest_tips.py
"""
from  tip_tip_distances import get_tip_tip_distances
from taxonomy_from_blast import get_common_lineages


class TaxonomyFromTipNeighbors(object):
    """get taxon for the nearest tips
    """
    def __init__(self, phylo_root):
        """
        - phylo_root: the root as PhyloNode of the phylo tree.
        """
        self._root = phylo_root

    def _init_tip_tip_distances(self):
        """ calc the tip_to_tip distance matrix
        """
        # ._tipidx added to each tip.
        self._dist, self._tip_order = get_tip_tip_distances(self._root)


    def _each_tip_neighbors(self, query_tips, data_tips,
            max_distance=None, band_width=None):
        """return [(tip, neighbors)].

        - query_tips: query
        - data_tips: tips to search
        - band_width=0: max delta distance with the best one.
        - max_distance=None: override band_width if not None, max distance from
          the query tip.
        """
        assert max_distance is None or band_width is None

        ## extract distance vectors
        row_idxs = [q._tipidx for q in query_tips]
        col_idxs = [d._tipidx for d in data_tips]
        dist = self._dist[row_idxs][:, col_idxs]

        ## calc the radiuses to draw the circle for neighbors
        if max_distance is not None:
            radiuses = repeat(max_distance, len(row_idxs))
        else:
            radiuses = dist.min(axis=1) + band_width

        ## draw each circle, yield the tip and neighbors
        data_tips = array(data_tips, object)
        for tip, row, radius in zip(query_tips, dist, radiuses):
            yield tip, data_tips[row <= radius]

    #wrapper
    def eachTipnameCommonLineage(self, tip_to_lineage, lineages_to_common,
            **tip_neighbor_kw):
        for tip, neighbors in self._each_tip_neighbors(tip_neighbor_kw):
            yield tip.Name, lineages_to_common(map(tip_to_lineage, neighbors))



## wrapper function for convenience
def get_taxonomy_from_tip_neighbors(phylo_tree, tipname_lineage_nds,
        taxonomy_tree, majority_threshold=1, **tip_neighbor_kw):
    """ return {tipname: common lineage}.

    - phylo_tree: phylogeny tree lines in NewWick format.
    - tipname_lineage_nds: Arb nds lines, each line as 'tipname\tlineage'
    #- taxonomy_tree: the taxonomy tree with valid names.
    - **tip_neighbor_kw: pass to ._each_tip_neighbor(,band_width=,
      max_distance=)
    """
    obj = TaxonomyFromTipNeighbors(DndParser(phylo_tree))
    tip_lineages = PhyloTaxTree._get_taxas_from_lines(tipname_lineage_nds)
    tip_to_lineage = lambda x: tip_lineages.get(x, [])

    query_tips = []; data_tips = []
    for tip in root.iterTips():
        if tip.Name in tip_lineages():
            data_tips.append(tip)
        else:
            query_tips.append(tip)

    lineages_to_common = partial(get_common_lineage,
            threshold=majority_threshold)
    return dict(obj.eachTipnameCommonLineage(tip_to_lineage,
        lineages_to_common, query_tips=query_tips, data_tips=data_tips,
        **tip_neighbor_kw)










