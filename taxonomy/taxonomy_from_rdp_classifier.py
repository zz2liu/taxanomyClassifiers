"""taxonomy_from_rdp_classifier.py

parse from rdp web assignment lines to each (seqname, taxon in domain pcofg
order)
"""
from itertools import dropwhile
from cogent.util.transform import first_index
from py_util.str_ import percentage
from util import RdpTaxonomy
from _rdp_taxonomy import RDP_TAXONOMY
from warnings import warn
warn('use parse_rdp.py instead.')

















def parse_tip_taxa(rdp_lines, rdp_taxonomy=None, threshold=0.0):
    """yield each (tipname, taxon as rdp titles)

    - rdp_lines: rdp classifier result txt lines
    - rdp_taxonomy=RDP_TAXONOMY: a util.RdpTaxonomy obj.
    - threshold=0.0: the min ratio of assignments to the tax to trust
      0.50 is a reasonable option.
    """
    threshold_percent = threshold * 100
    get_titles = RdpTitles(rdp_taxonomy or RDP_TAXONOMY)
    
    ## throw away the header lines
    assign_lines = dropwhile(lambda x: x != 'Details:\n', rdp_lines)
    assign_lines.next() #ignore the 'Details:' line
    
    for line in assign_lines:
        tipname, taxon, percentages = _get_tip_lineage(line)
        first_badidx = first_index(percentages, lambda x: x < threshold_percent)
        yield tipname, get_titles(taxon[:first_badidx])
        
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
    
def RdpTitles(rdp_taxonomy):
    """return a fun(lineage) -> rdp_titles"""
    def fun(lineage):
        tax_array = rdp_taxonomy.taxArrayFromTaxon(lineage)
        result = rdp_taxonomy.taxonFromTaxArray(tax_array)
        return result
    return fun
        
