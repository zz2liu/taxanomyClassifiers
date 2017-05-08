"""ncbi_gi_to_lineage.py 

06/25/07 created offline version first.
todo: make a dmp like that of tip_taxa.

"""
import shelve
from cogent.parse.ncbi_taxonomy import TaxonomyFromFiles

def _make_shelf(dmp_path, out_shelf_path):
    shelf = shelve.open(out_shelf_path, 'c')
    for line in open(dmp_path):
        gi, taxid = line.split()
        shelf[gi] = taxid
    shelf.close()

def each_gi_to_taxid(gis, shelf_path):
    shelf = shelve.open(shelf_path, 'r')
    for gi in gis:
        yield shelf[gi]

def each_taxid_to_lineage(taxids, nodes_dmp, names_dmp):
    """yield each lineage as [(rank, name),] in leaf-root order.

    Note: both rank and name can be duplicated, some ranks are just 'no rank'.
    """
    tx = TaxonomyFromFiles(file(nodes_dmp), file(names_dmp))
    for taxid in taxids:
        node = tx.ById[taxid]
        yield [(p.Rank, p.Name) for p in n.ancestors()]

def each_gi_to_titles(gis, ranks=['phylum', 'class', 'order', 'family',
'genus'], default_name='Unknown'):
    """hard coded for now.

    by 'titles', I mean [phylum_name, class_name, ...]
    """
    shelf_path = ''
    nodes_dmp = ''
    names_dmp = ''

    each_taxid = each_gi_to_taxid(gis, shelf_path)
    each_lineage = each_taxid_to_lineage(each_taxid, nodes_dmp, names_dmp)
    for lineage in each_lineage:
        rank_names = dict(lineage)
        titles = [rank_names.get(rk, default_name)
                for rk in ranks)
        yield titles


    

