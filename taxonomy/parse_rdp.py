"""taxonomy_from_rdp_classifier.py

parse the rdp classification results to a list of (d,p,c,o,f,g)
"""
from warnings import warn
from itertools import dropwhile
from time import time
from pdb import set_trace
from pickle import dump, load
from cogent.util.transform import first_index

from py_util.str_ import percentage
from py_util.tree_new import DndParser
from taxonomy.util import Taxonomy
from taxonomy.data import BERGY78_TAX_TAB
from taxonomy.lineage_tree import make_getter

def each_lineage(rdp_lines, min_percentage=0):
    """yield each (tipname, taxon as rdp titles in domain, ..., genus order)

    Note: the dpcofg list is of variable lengths.
    - rdp_lines: rdp classifier result txt lines
    - dpcofg_getter: a func(lineage as list) -> dpcofg list
    - threshold=0.0: the min ratio of assignments to the tax to trust
      0.50 is a reasonable option.
    """
    ## throw away the header lines
    assign_lines = dropwhile(lambda x: not x.startswith('Details:'), rdp_lines)
    assign_lines.next() #ignore the 'Details:' line
    
    for line in assign_lines:
        if not line.strip(): continue
        tipname, taxon, percentages = _get_tip_lineage(line)
        first_badidx = first_index(percentages, lambda x: x < min_percentage)
        yield tipname, taxon[:first_badidx]
        
def _get_tip_lineage(line):
    """return (tipname, lineage as str list, percentages as a list),
    starting from domain level (hard coded !).

    - line: a rdp classifier assignment line. like
    '16254;  ; Root; 100%; Bacteria; 100%; Cyanobacteria; 89%;'\
        'Cyanobacteria; 89%; Family 1.1; 54%; Prochlorococcus; 24%'
    """
    fields = line.strip().split('; ')
    tipname = fields[0]
    lineage = fields[4::2]
    percentages = map(percentage, fields[5::2])
    return tipname, lineage, percentages
    
       
# the main horse
def parse(rdp_lines, tax_lines=None, using_dpcofg=False, min_percentage=50):
    """return a dict of {seqname: [dpcofg]}

    - rdp_lines: the ori result file from rdp classifier
    - tax_lines: the ori bergy tax tab file
    - using_dpcofg=False: two methods to get dpcofg from ori lineage.
    """
    if not tax_lines:
        tax_lines = open(BERGY78_TAX_TAB)
    if using_dpcofg:
        #build a Taxonomy tree -- to be improved (unnecessary complicated)
        genus_tree = build_genus_tree(tax_lines)
        dpcofg_tree = RdpTaxonomy.fromRawTaxTree(genus_tree)
        # __str__ of TreeNode is not reversible
        getter = FromDpcofgTree(dpcofg_tree._root.str().splitlines())
    else: #useing ori tax tab
        getter = make_getter(tax_lines)

    return list((seqname, getter(lineage)) for seqname, lineage in
            each_lineage(rdp_lines, min_percentage=min_percentage))

def parse_many(in_fnames, out_fnames=None, taxlines=None,
        out_suffix='.dpcofg_list', min_percentage=50):
    """parse many rdp classification result files, return the times used.
    """
    out_fnames = out_fnames or [name + out_suffix for name in in_fnames]
    result = []
    for in_name, out_name in zip(in_fnames, out_fnames):
        print '%s -> %s ...' % (in_name, out_name)
        start = time()
        res = parse(open(in_name), taxlines, min_percentage=min_percentage)
        dump(res, open(out_name, 'w'))
        result.append(time()-start)
    return result

def main(argv):
    """usage example:

    python parse_rdp.py classifictations.txt --min_percentage=50,  > tmp.list
    """
    from py_util.cmdline import run, wrap
    run(argv, {
        None: wrap([open, open],
            {'min_percentage': float}, dump
            )(parse),
        'many': wrap(['safe_eval', 'safe_eval'],
            {'min_percentage': float})(parse_many),
        })








#deprecated
def FromDpcofgTree(dpcofg_tre_file):
    """return a fun(lineage) -> rdp_titles,
    
    - taxonomy_file: a file/lines of a strict dpcofg tree.
    """
    taxonomy_root = DndParser(dpcofg_tree.str(), Taxonomy)

    def fun(lineage):
        tax_node = taxonomy_root.taxNodeFromTaxon(lineage, False) #w/o root
        if tax_node is None:
            raise ValueError('Not found in taxonomy tree for lineage: %s'
                    % lineage)
        return tax_node.getTaxon()[1:] #w/o root
    return fun
if __name__ == '__main__':
    import sys
    main(sys.argv)
