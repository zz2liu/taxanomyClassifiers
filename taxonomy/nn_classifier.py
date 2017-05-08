"""nn_classifier.py  - Nearest Neighbor Classifier
"""

class NNClassifier(object):
    def __init__(self, dist_f, distmat_f=None):
        self._dist_f = dist_f
        self._distmat_f = distmat_f

    def train(self, std_ids, std_lineages, std_rows, std_distmat=None):
        """calc .std_distmat."""

    def classifyOne(self, query_row, neigbor_thres=None, max_diff=0,
            min_maj=1):
        """return a predicted lineage from nearest neibors."""
        dist_f = self._dist_f
        distances = [dist_f(query_row, std_row)
                for std_row in std_rows]
        #nn_idxs = nn_f(distances)
        nn_lineages = std_lineages[nn_idxs]
        #result_lineage = lineage_reduce_f(nn_lineages)
        return result_lineage

    def classifyMany(self, query_ids, query_rows):
        """return [(query_id, result_lineage),].
        """



