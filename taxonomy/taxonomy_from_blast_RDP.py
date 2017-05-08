"""taxonomy_from_blast_RDP.py
"""
from pdb import set_trace
from os import path

from py_util import Table
from taxonomy_from_tipnames import PhyloTaxTree, RdpTaxonomy
from taxonomy_from_blast import (each_query_besthits, get_common_lineage,
        SUBJECT_ID)


TAX_DIR = '/home/zongzhi/Projects/misc/zongzhi/taxonomy'
CORESET_RDP_NDS_PATH = path.join(TAX_DIR, 'coreset_rdp_good.nds')
get_first_word = lambda line: line.split()[0]

### from tipnames
def coreset_id_to_rdp_taxa(coreset_rdp_nds):
    """return {coreset_id: [rdp tax] in order of levels}
    
    - coreset_rdp_nds: lines of an Arb nds file, each line as
      'coreset_id\tRDP lineage'

    Warn: hard coded for now.
    """
    rdp_tax_tree_path = path.join(TAX_DIR, 'RDP_taxonomy.tre')
    rdp_taxonomy = RdpTaxonomy(file(rdp_tax_tree_path))

    coreset_id_taxon_table = Table(coreset_rdp_nds, header=False)
    result = {} #{coreset_id: [rdp_tax]}
    for name, taxon in coreset_id_taxon_table:
        if not taxon: continue  #some id with no taxon defined
        tax_array = rdp_taxonomy.taxArrayFromTaxon(taxon.split('; '))
        result[name] = rdp_taxonomy.taxonFromTaxArray(tax_array)
    return result

### from blast
def each_query_top_coresets(in_blast):
    """return {query: top_coresets as list}
    """
    coreset_id_from_hit = lambda line: line.split('greengenes|')[1]
    hit_parser=lambda row: coreset_id_from_hit(row[SUBJECT_ID])
    return each_query_besthits(in_blast, hit_parser=hit_parser,
            query_parser=get_first_word)


### a wrapper combining from_blast and from_tipnames
def get_query_common_lineages(in_blast, coreset_rdptaxa):
    """return {query: common_lineage}

    - in_blast: lines of blast9 result.
    - coreset_rdp_taxa: {coreset_id: rdp_taxa}
      or lines of the Ard nds file. 'coresetname\trdp_taxon'
    """
    if not isinstance(coreset_rdptaxa, dict):
        coreset_rdptaxa = coreset_id_to_rdp_taxa(coreset_rdptaxa)
    result = {}
    for query, top_coresets in each_query_top_coresets(in_blast):
        top_lineages = [coreset_rdptaxa.get(id, [])
                for id in top_coresets]
        result[query] = get_common_lineage(top_lineages)
    return result

if __name__ == "__main__":
    ## quick test
    res = coreset_id_to_rdp_ranks()
    set_trace()

 



