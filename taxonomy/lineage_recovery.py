"""










Deprecated?
tobe rewritten using working.taxononmy_finishing modules
"""
from __future__ import division
from numpy import *
from py_util.dict_ import calc_freqs
from warnings import warn


def get_lineage_recovery(obs, exp):
    """return a vector of 1/0/-1, representing tax recovered or not or na.

    - obs, exp: obs and exp lineage vector of tax or None.
      they should be the same length (no check).
    """
    warn('deprecated, use')
    result = []
    for o, e in zip(obs, exp):
        if o and e:
            result.append(o == e)
        else:
            result.append(-1)
    return result

def calc_recovery_scores(recoveries):
    """return tax recovery proportions and covery proportions.

    - recoveries: a seq by level 2dlist of 1/0/-1, with 1: recoved, 0: not
      recovered and -1: NA.
    """
    recoveries = array(recoveries, int)
    result = []
    for col in recoveries.T:
        valid_col = col[col >= 0]
        accuracy = valid_col.sum()/len(valid_col)
        coverage = len(valid_col)/len(col)
        result.append((accuracy, coverage))
    return zip(*result)




######
# not tested
def get_lineage_recoveries(obs_lineages, exp_lineages):
    """yield each bool list of tax recoveries."""
    for (obs, exp) in zip(obs_lineages, exp_lineages):
        yield get_lineage_recovery(obs, exp)

def _get_query_common_lineages(query_hitids, dpcofg, threshold=1,):
    """return {query: common_lineage}

    - dpcofg: {seq_id: [domain,phylum, c, o, f, g]}
    """
    result = {}
    for query, top_hits in query_hitids:
        top_lineages = [dpcofg[id] for id in top_hits]
        result[query] = get_common_lineage(top_lineages, threshold=threshold)
    return result

def get_query_common_lineages(in_blast, dpcofg, max_rel_diff=0, threshold=1,
        exclude_single_generus_query=False):
    """return {query: common_lineage}

    - in_blast: lines of blast9 result.
    - dpcofg: {seq_id: [domain,phylum, c, o, f, g]}
    """
    query_hitids = query_besthits_without_self(in_blast,
            max_rel_diff=max_rel_diff)

    if exclude_single_generus_query:
        single_ids = get_single_genus_ids(dpcofg)
        for i in single_ids:
            del query_hitids[i]

    result = _get_query_common_lineages(query_hitids, dpcofg,
            threshold=threshold)
    return result

def get_query_recoveries(query_lineages, expect_lineages):
    """return a bool 2dlist of tax recoveries and query names.

    inputs exp to be dicts.
    """
    query_ids = query_lineages.keys()
    result = []
    for k in query_ids:
        obs, exp = query_lineages[k], expect_lineages[k]
        lacks = len(exp) - len(obs)
        obs = obs + [None]*lacks
        result.append(get_lineage_recovery(obs, exp))
    return result, query_ids

def get_single_genus_ids(dpcofg):
    """dpcofg: a {id: [domain, phylum, ...]}
    """
    genera = [v[-1] for v in dpcofg.values()]
    genus_freqs = calc_freqs(genera)
    result = [id for id, lin in dpcofg.iteritems()
            if genus_freqs[lin[-1]] == 1]
    return result
