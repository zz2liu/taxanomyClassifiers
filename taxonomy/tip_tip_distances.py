"""tip_tip_distances.py - using dynamic programming!

Functions to give each tip_to_tip_distance for a tree.

07/16/07 tip_order change to [tipnode]
08/13/07 assign_tip_range seperated from the original function.
02/27/08 use add.outer to speedup.

todo: make .assignTipRange to TreeNode?
todo: mv to tree_new.py?
"""
from py_util.cbook import comb, cross_comb
from numpy import zeros, add
from itertools import count
from py_util.util import disp, set_trace
add_outer = add.outer

def assign_tip_range(root, attr_name='_tip_range'):
    """assign tip_range to root and each_descendents.

    - root: a TreeNode obj as root of a tree.
    - attr_name: the attr name to be set for tip_range.
    """
    for i, tip in enumerate(root.iterTips(mode='postorder')):
        setattr(tip, attr_name, (i, i+1))

    for node in root.iterNontips(mode='postorder'):
        tip_start = getattr(node.Children[0], attr_name)[0]
        tip_stop = getattr(node.Children[-1], attr_name)[1]
        setattr(node, attr_name, (tip_start, tip_stop))

def get_tip_tip_distances(root, verbose=False):
    """return a distance matrix between tips and a tip order.

    - root: a TreeNode as root of the tree.
    Note: ._tiprange added to each node.
    """
    ## linearize the tips in postorder.
    # ._start, ._stop compose the slice in tip_order.
    tip_order = list(root.iterTips())
    num_tips = len(tip_order)
    result = zeros((num_tips, num_tips), float) #tip by tip matrix
    tipdistances = zeros((num_tips), float) #distances from tip to curr node

    def update_result(): #using i node, i tipdistances, o result
        """set tip_tip distance between tips of different child"""
        for child1, child2 in comb(node.Children, 2):
            s1, s2 = slice(*child1._tip_range), slice(*child2._tip_range)
            curr = add_outer(tipdistances[s1], tipdistances[s2])
            result[s1, s2] = curr
            result[s2, s1] = curr.T

    assign_tip_range(root, '_tip_range')
    node_counter = count()
    for node in root.iterNontips(mode='postorder'):
        if verbose: disp('\r%i   ' % node_counter.next())
        for child in node.Children:
            s = slice(*child._tip_range)
            tipdistances[s] += child.Length
        ## update result if possible
        if len(node.Children) > 1: #not single child
            update_result()
    return result, tip_order











#### obsoleted functions below
def get_tip_tip_distances_dict(root):
    """return {(tip, other_tip): distance}.

    - root: a TreeNode as root of the tree.
    """
    ## tipdistances: {tipname: distance to the curr node}
    result = {} #{(tip, other_tip): distance}
    for node in root.iterNodes(mode='postorder'):
        ## init miminal subtree
        if node.isTip():
            node._tipdistances = {node.Name: 0}
            continue
        ## subtree with solved child wedges
        node._tipdistances = {} #tobe node.tipdistances
        for child in node.Children:
            for tip, distance in child._tipdistances.iteritems():
                node._tipdistances[tip] = distance + child.Length

        ## update result if nesseccery
        if len(node.Children) == 1: continue #single child
        for child1, child2 in comb(node.Children, 2):
            child1_tips = list(child1._tipdistances.keys())
            child2_tips = list(child2._tipdistances.keys())
            for tip1, tip2 in cross_comb([child1_tips, child2_tips]):
                result[(tip1, tip2)] = \
                        node._tipdistances[tip1] +\
                        node._tipdistances[tip2]
    return result




def get_tip_tip_distances_array(root):
    """return a distance matrix between tips and a tip order.

    - root: a TreeNode as root of the tree.

    Warning: ._start and ._stop added to root and its descendents.
    """
    ## linearize the tips in postorder.
    # ._start, ._stop compose the slice in tip_order.
    tip_order = list(root.iterTips(mode='postorder'))
    for i, tip in enumerate(tip_order):
        tip._start, tip._stop = i, i+1

    num_tips = len(tip_order)
    result = zeros((num_tips, num_tips), float) #tip by tip matrix
    tipdistances = zeros((num_tips), float) #distances from tip to curr node

    def update_result(): # set tip_tip distance between tips of different child
        for child1, child2 in comb(node.Children, 2):
            for tip1 in range(child1._start, child1._stop):
                for tip2 in range(child2._start, child2._stop):
                    result[tip1, tip2] = tipdistances[[tip1, tip2]].sum()

    for node in root.iterNontips(mode='postorder'):
        ## subtree with solved child wedges
        starts, stops = [], [] #to calc ._start and ._stop for curr node
        for child in node.Children:
            tipdistances[child._start : child._stop] += child.Length
            starts.append(child._start); stops.append(child._stop)
        node._start, node._stop = min(starts), max(stops)
        ## update result if nessessary
        if len(node.Children) > 1: #not single child
            update_result()
    return result+result.T, tip_order

def get_tip_tip_distances_array_old(root):
    """return a distance matrix between tips and a tip order.

    - root: a TreeNode as root of the tree.

    Warning: ._tip_idx and ._tips add to root and its descendents.
    """
    tip_order = list(root.iterTips(mode='postorder'))
    for i, tip in enumerate(tip_order):
        tip._tip_idx = i

    num_tips = len(tip_order)
    result = zeros((num_tips, num_tips), float) #tip by tip matrix
    tipdistances = zeros((num_tips), float) #distances from tip to curr node

    def update_result(): # set tip_tip distance between tips of different child
        for child1, child2 in comb(node.Children, 2):
            for tip1, tip2 in cross_comb((child1._tips, child2._tips)):
                result[tip1, tip2] = tipdistances[tip1] + tipdistances[tip2]

    for node in root.iterNodes(mode='postorder'):
        ## init miminal subtree (tip)
        if node.isTip():
            node._tips = [node._tip_idx]
            continue
        ## subtree with solved child wedges
        node._tips = [] #tips of this node as indexes
        for child in node.Children:
            tipdistances[child._tips] += child.Length
            node._tips.extend(child._tips)
        ## update result if nessessary
        if len(node.Children) > 1: #not single child
            update_result()
        ## free the _tips no more useful for space performance.
        for child in node.Children:
            del child._tips
    return result+result.T, tip_order

