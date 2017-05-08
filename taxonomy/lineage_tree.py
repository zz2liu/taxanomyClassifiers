"""lineage_tree.py - Provide functions to convert lineage str to dpcofg.

Similar Code:
    taxonomy_tree.py
    taxonomy_from_rdp_classifier.py
    my_util.tree_new.search

03/04/08 pad_mid_missing_tax improved

todo: clean the code, make orthoganal functions, reduce the deps.
"""
from py_util.parse.tree import tree_from_tab
from py_util.tree_new import TreeNode
from pdb import set_trace

RANKS = {'domain':0, 'phylum':1, 'class':2, 'order':3, 'subsection':3, #rel5.1
        'family':4, 'genus':5, 'species': 6, 'subclass': 2.5, 'suborder': 3.5,
        'no rank': -1, 'norank': -1}

def get_rank_number(rank_name):
   rank_name=rank_name.strip()
   return RANKS[rank_name]

def build_genus_tree(tax_tab, **kw):
    """make tax_tree deep to genus rank with tax_tab_lines.
    **kw: pass to tree_from_tab(, id_idx, name_idx, parent_idx, root_id=,
    sep=, ...)
    """
    kw.setdefault('other_idxs', {'rank':-1})
    kw.setdefault('sep', '*')
    kw.setdefault('field_constructors',
            [int, str, int, int, get_rank_number])
    kw.setdefault('ignore', lambda x: x[-1] == 6)
    return tree_from_tab(tax_tab, 0,1,2, **kw)

def lineage_to_dpcofg(tree, lineage, num_ranks=6):
    """return a list of tax/'' in order of dpcofg.

    - tree: root node of the tax tree with rank numbers.
      using .Name, .params['rank'].
    - lineage: the list of tax names without root name.

    Note: the result lineage will not be garanteed to get the same result if
    you input to this function again, especially when there is gap or subclass
    in the original lineage.
    """
    ## walking the tree with lineage, collecting {.rank: .name}
    result = {}
    wedge = tree
    for name in lineage:
        for n in wedge.Children:
            if n.Name == name:
                result[n.params['rank']] = name
                wedge = n
                break #continue to the next name
        else: #not found in wedge's children
            break
    ## return a list of dpcofg out of the {rank: name}
    return [result.get(i, '') for i in range(num_ranks)]

def make_getter(tax_tab_file, using_dpcofg=False, **kw):
    """return a func(lineage) -> dpcofg with unknown as ''.

    a quick wrapper of lineage_to_dpcofg.
    _
    - tax_tab_lines: a tab sep lines with ranks.
    tax_tree = build_genus_tree(tax_tab_lines)

    fvs: build_genus_tree, lineage_to_dpcofg
    """
    if isinstance(tax_tab_file, TreeNode):
        genus_tree = tax_tab_file
    else:
        genus_tree = build_genus_tree(tax_tab_file, **kw)

    if using_dpcofg:
        raise NotImplementedError()
    else:
        def getter(lineage):
            """return dpcofg from lineage start from domain."""
            result = lineage_to_dpcofg(genus_tree, lineage)
            return [tax or '' for tax in result]
    return getter

##no dedicated test yet
#def pad_mid_missing_tax_old(dpcofg):

    #idxs = [] # to be pad
    #flag = False
    #for i in range(6)[::-1]:
        #if dpcofg[i]:
            #flag = True
            #continue
        #if flag and not dpcofg[i]:
            #idxs.append(i)
    #result = dpcofg[:] #copy
    #for i in idxs[::-1]:
        #result[i] = result[i-1] + '_'
    #return result

def pad_mid_missing_tax(dpcofg):
    """pad the missing_tax in the middle with ancestor name + '_'

    - dpcofg: a list of tax names in the order of dpcofg, unknow tax as ''.
    """
    for i in range(len(dpcofg))[::-1]:
        if dpcofg[i]:
             break
    for j in range(1, i):
        if not dpcofg[j]:
            dpcofg[j] = dpcofg[j-1] + '_'
    return dpcofg

def each_id_lineage(id_lineage_file, field_sep='\t',lineage_sep=';',
        root_name='root'):
    """yield eadh (seqname, lineage as raw list of taxnames"""
    root_name = root_name.lower()
    for line in id_lineage_file:
        name, lineage = line.strip().split(field_sep)
        lineage = map(str.strip, lineage.split(lineage_sep)) #rip off root
        if lineage and lineage[0].lower() == root_name:
            lineage = lineage[1:]
        yield name, lineage

def each_seqname_dpcofg(tax_tab_file, id_lineage_file, **kw):
    """yield each (seqname, dpcofg as 6-item list).

    - tax_tab_file: eg, open(data.BERGY78_TAX_TAB)
    - id_lineage_file: eg, open(data.BERGY78_LINEAGES)
    **kw:
     field_sep='\t', lineage_sep=';', root_name='root'

    fvs: make_getter, pad_mid_missing_tax
    """
    if callable(tax_tab_file):
        get_dpcofg = tax_tab_file
    else:
        get_dpcofg = make_getter(tax_tab_file)

    for name, lineage in each_id_lineage(id_lineage_file, **kw):
        dpcofg = get_dpcofg(lineage)
        # must pad before get common lineage
        yield name, pad_mid_missing_tax(dpcofg)

def build_dpcofg_dict(tax_tab_file, id_lineage_file, **kw):
    """return a dict of id -> dpcofg as list.
    """
    return dict(each_seqname_dpcofg(tax_tab_file, id_lineage_file, **kw))
build_dpcofg_dict.__doc__ += '\n'.join(
        each_seqname_dpcofg.__doc__.splitlines()[2:])








#deprecated, use make_getter instead, for a better name though

