#!/usr/bin/env python
"""
>>> from cogent.parse.ncbi_taxonomy import TaxonomyFromFiles
>>> tx = TaxonomyFromFiles(file('nodes.dmp'), file('names.dmp'))
>>> n = tx.ById[100]
>>> [(p.Name, p.Rank) for p in n.ancestors()]
[('Ancylobacter', 'genus'), ('Xanthobacteraceae', 'family'), ('Rhizobiales', 'order'), ('Alphaproteobacteria', 'class'), ('Proteobacteria', 'phylum'), ('Bacteria', 'superkingdom'), ('cellular organisms', 'no rank'), ('root', 'no rank')]
>>> n = tx.ByName['Proteobacteria']
>>> [(p.Name, p.Rank) for p in n.ancestors()]
[('Bacteria', 'superkingdom'), ('cellular organisms', 'no rank'), ('root', 'no rank')]

WorkFlow: 
    fasta files
    blast against nt (fasta_file) -> {fasta_id: [hit_ids]}
    fastacmd -T(unique_hit_ids) -> {hit_id: tax_id}
    NcbiTaxonomy(unique_tax_ids) -> {tax_id: lineage}
    ==
    -> {fasta_id: [lineages]}

    -> {fasta_id: taxon}

Todo: quality control in each step?
Todo: use the decorator autopickle

"""
from __future__ import division
from functools import partial
from pdb import set_trace
from os import popen, path, system
from glob import glob
from itertools import chain
from cPickle import load, dump
from py_util.db.ncbi import NcbiTaxonomy
from numpy import concatenate, array

from cogent.parse.blast import MinimalBlastParser9, FastacmdTaxonomyParser
from cogent.db.ncbi import taxon_ids_to_lineages
from cogent.maths.stats.util import UnsafeFreqs as Freqs
from cogent.parse.fasta import MinimalFastaParser

from py_util.table import Table
from py_util.str_ import extract_between
from py_util.file_ import File
from py_util.misc import lace #map(None, *lists) fails with single item lists.

SUBJECT_ID = 1
BIT_SCORE = -1
WORKING_DIR = '.'
DEFAULT = True

####
# get_common_lineage and support funs.
def get_common_lineage(lineages, threshold=1, filter_=None):
    """return the common lineage in root first order.

    - lineages: a list of lineages, each as a title list in the same
        order.
    - threshold: how similar can we say 'common'.  (0.5, 1]

    deprecated(filter_)
    decision: keep it simple.
    todo: divide the second part as get_dominant_lineage?
    """
    if not 0.5 <= threshold <= 1:
        raise ValueError('expected 0.5 <= threshold <= 1, got %s'
            % threshold)

    if filter_:
        lineages = filter(filter_, lineages)

    if not lineages:
        return [None]
    if len(lineages) == 1: #single lineage ['bact', '', '', 'unknown']
        return lineages[0]
    if threshold == 1:
        return _get_common_heads(lineages)
    else: #if 0.5<threshold<1:
        return _get_dominant_heads(lineages, threshold)

def _get_dominant_heads(lineages, threshold=0):
    result = []
    for taxs in lace(*lineages):
        freqs = Freqs.newFromSeq(taxs)
        mode = freqs.Mode
        if mode is None or freqs[mode]/freqs.Sum < threshold:
            break
        result.append(mode)
    return result

def _get_common_heads(lineages):
    result = []
    for taxs in lace(*lineages): #map(None, *lineages):
        if len(set(taxs)) == 1:# and taxs[0]: #common this level
            result.append(taxs[0])
        else: #different at this level
            break
    return result


#####
# each_query_besthits and extensons
def _each_query_besthits(in_blast, hit_parser=DEFAULT, query_parser=None):
    """return each (query, [lineages]).

    - in_blast: blast9 result as lines.
    - hit_parser: fun(row) -> hit_id
    """
    if hit_parser is DEFAULT:
        hit_parser = lambda row: row[SUBJECT_ID]
    for props, data in MinimalBlastParser9(in_blast):
        query = props['QUERY']
        if query_parser: query = query_parser(query)
        
        top_score = data[0][BIT_SCORE]
        top_hits =  [hit_parser(row) for row in data
                if row[BIT_SCORE] == top_score]
        yield query, top_hits
each_query_besthits = _each_query_besthits #alias for export

gi_from_line = lambda line: int(extract_between(line, 'gi|', '|'))
each_query_best_gis = partial(_each_query_besthits,
        hit_parser=lambda row: gi_from_line(row[SUBJECT_ID]))


def each_query_toplineages(in_blast, seqid_lineages):
    for query, top_hits in _each_query_besthits(in_blast):
        top_lineages = [seqid_lineages[hit] for hit in top_hits]
        yield query, top_lineages



#class not tested yet
class TaxonomyFromBlast(object):
    """A work f low class, works as a function organizer.

    todo: generalize dump methods.
    """
    def __init__(self):
        self.BlastLines = None
        self.QueryGis = self.GiTaxids = self.TaxidLineages = None

    def getQueryGis(self):
        """return and dump .QueryGis - {query_id: [hit_gis]}
        """
        dump_path = '%s/BlastToLineages.QueryGis.dump' % WORKING_DIR
        if not self.BlastLines:
            raise DependError('.BlastLines not given')

        try:
            result = load(file(dump_path))
        except IOError:
            result = dict(each_query_gis(self.BlastLines))
            dump(result, file(dump_path, 'w'))

        self.QueryGis = result
        return result

    def getGiTaxids(self):
        """return and dump .GiTaxids - {unique_GI: taxid}
        """
        dump_path = '%s/BlastToLineages.GiTaxids.dump' % WORKING_DIR
        if not self.QueryGis:
            self.QueryGis = self.getQueryGis()

        try:
            result = load(file(dump_path))
        except IOError:
            unique_gis = set(concatenate(self.QueryGis.values()))
            result = dict(each_gi_taxonId(unique_gis))
            dump(result, file(dump_path, 'w'))

        self.GiTaxids = result
        return result

    def getTaxidLineages(self):
        """return and dump .TaxidLineages - {unique_taxid: lineage}
        """
        dump_path = '%s/BlastToLineages.TaxidLineages.dump' % WORKING_DIR
        if not self.GiTaxids:
            self.GiTaxids = self.getGiTaxids()

        try:
            result = load(file(dump_path))
        except IOError:
            unique_taxids = set(self.GiTaxids.values())
            result = dict(each_taxid_lineage(unique_taxids))
            dump(result, file(dump_path, 'w'))

        self.TaxidLineages = result
        return result

    def eachQueryLineages(self):
        """work with the three dicts.
        """
        query_gis = self.getQueryGis()
        gi_taxids = self.getGiTaxids()
        taxid_lineages = self.getTaxidLineages()

        for query, gis in query_gis.items():
            taxids = map(gi_taxids.__getitem__, gis)
            lineages = map(taxid_lineages.__getitem__, taxids)
            yield query, lineages

    @classmethod
    def eachQueryCommonTaxon(cls, query_lineages, **kw):
        """return a dict of {query_id: [taxs]}
        """
        for query, lineages in query_lineages:
            common_lin = cls.get_common_lineage(lineages, **kw)
            #print common_lin
            yield query, common_lin

    get_common_lineage = staticmethod(get_common_lineage)

    




def parse_gi(line):
    """return the gi  str from ncbi seq id line

    example:
    'gi|89112992|gb|DQ401104.1|' -> '89112992'
    """
    return line.split('|')[1]

def each_query_gis(blast_lines, hit_parser=parse_gi):
    """yield each que ry_id and the set of hit GIs from blast.out.

    todo: remove frozenset, keep the original order and scores?
    """
    for props, data in MinimalBlastParser9(blast_lines):
        yield props['QUERY'], frozenset(hit_parser(row[SUBJECT_ID]) for row in data)

def each_gi_taxonId(seq_gis):
    """yield each (gi,  taxon_id) for gis provided.

    - seq_gis: gi numbers as lines.

    Warning: fastacmd line hard coded for now.
    """
    #tmp = as_filename('\n'.join(map(str, seq_gis)))
    tmp = File.fromLines(seq_gis, close=True).name
    result = popen('fastacmd -d env_nt -i %s -T' % tmp)
    for rec in FastacmdTaxonomyParser(result):
        yield parse_gi(rec['seq_id']), rec['tax_id']

def each_taxid_lineage(tax_ids):
    """return each (taxid, lineage).
    """
    for taxid, name, lineage in NcbiTaxonomy().eachIdNameLineage(tax_ids):
        yield taxid, lineage

#### quick codes
def each_seqid_lineage_from_european_rdp(embl_files):
    """return each (seqid, lineage as a list) from embl files.

    - embl_files: a list of file objects in embl format.
    """
    from snp.ebi_like import Parser, labeloff

    def to_str(lines):
        return ' '.join([line[5:] for line in lines])

    for rec in Parser(chain(*embl_files)):
        id = to_str(rec['ID'])
        oc = to_str(rec['OC'])
        seqid = '_'.join(id.split())
        taxon = map(str.strip, oc.split(';'))
        yield seqid, taxon

def each_query_lineages(blast_files, seqid_lineages):
    """return each (query, [lineages]).

    - seqid_lineages: {seqid: lineage}
    """
    for props, data in MinimalBlastParser9(chain(*blast_files)):
        query = props['QUERY']
        hits =  [row[SUBJECT_ID] for row in data]
        lineages = []
        for h in hits:
            lineages.append(seqid_lineages[h])
        yield query, lineages


#####
# filepaths
working_dir = '/home/zongzhi/Working/noah'
seqid_lineages_dump = './data/seqid_lineages_european_rrna.dump'

embl_files = ['%s/%s.embl' % (working_dir, fn)
            for fn in ['arch', 'bact', 'euka', 'mito', 'plastid']]
fastafiles = glob('%s/fasta/*.fasta' % working_dir)
blastfiles = [('%s.blast' % n) for n in fastafiles]
query_lineages_dump = '%s/query_lineages.dump' % working_dir
query_toplineages_dump = '%s/query_toplineages.dump' % working_dir
query_lineage_table = '%s/query_common_top_1.xls' % working_dir




def nearest_common_lineages_from_blast(in_blast, seqid_lineages,
        **kw):
    """workflow: return [(seqid, common_lineage)].

    **kw: pass to .get_common_lineage(,threshold=1,)
    """
    result = []
    for query, top_lineages in each_query_toplineages(
            in_blast, seqid_lineages):
        id = query.split()[0]
        common_lineage = TaxonomyFromBlast.get_common_lineage(
                top_lineages, **kw)
        result.append((id, common_lineage))
    return result


#to rm
###adding midwest

########run codes
def run_dump_seqid_lineages():
    files = map(file, embl_files)
    result = list(each_seqid_lineage_from_european_rdp(files))
    dump(result, file(seqid_lineages_dump, 'w'))

def run_dump_query_lineages():
    seqid_lineages = dict(load(file(seqid_lineages_dump)))
    files = map(file, blastfiles)
    result = list(each_query_lineages(files, seqid_lineages))
    dump(result, file(query_lineages_dump, 'w'))

def run_dump_query_toplineages():
    seqid_lineages = dict(load(file(seqid_lineages_dump)))
    files = map(file, blastfiles)
    result = list(each_query_toplineages(files, seqid_lineages))
    dump(result, file(query_toplineages_dump, 'w'))

def run_query_common_taxon(query_lineages,  out_xls, threshold=1,):
    #query_lineages  = load(file(query_lineages_dump))
    tb = TaxonomyFromBlast()
    result = list(tb.eachQueryCommonTaxon(query_lineages, threshold=threshold,
            filter_=None))#lambda x: x and len(x)>4))
    #write to xls
    outfile = file(out_xls, 'w')
    for query, taxon in result:
        if not taxon: taxon = []
        outfile.write('%s\t%s\n' % (query, '\t'.join(taxon)))
    outfile.close()

def run_blast_all(fastafiles, blastfiles):
    for infile_name, outfile_name in zip(fastafiles, blastfiles):
        system('blastall -p blastn -e 1e-10 -b 10 -m 9 -d all_ssu '
                '-i %s -o %s' % (infile_name, outfile_name))

def run_dump_nearest_common_lineages(in_fasta, in_blast, out_xls=None):
    """using seqid_lineages_dump,"""
    query_clones = get_query_clones_from_midwest_fasta(in_fasta)
    print query_clones.items(), len(query_clones)

    seqid_lineages = dict(load(file(seqid_lineages_dump)))
    query_lins = dict(nearest_common_lineages_from_blast(
        in_blast, seqid_lineages))
    print query_lins.items()[:5], len(query_lins)

    clone_lins = dict((clone, query_lins[q])
            for q, clone in query_clones.items())
    print clone_lins.items()[:5], len(clone_lins)
    dump(clone_lins, file(midwest_clone_lineages_dump, 'w'))

def run_make_bact_fungi_other_fasta(in_fasta, out_bact, out_fungi, out_other):
    clone_lineages = load(file(midwest_clone_lineages_dump))
    for label, seq in MinimalFastaParser(in_fasta):
        query, desc = label.split(None, 1)
        clone = desc.split()[3]
        lin = clone_lineages[clone]
        fasta = '>%s %s\n%s\n' % (clone, '; '.join(lin), seq)
        if lin[0:1] == ['Bacteria']:
            out_bact.write(fasta)
        elif lin[1:2] == ['Fungi']:
            out_fungi.write(fasta)
        else:
            print lin
            out_other.write(fasta)
    out_bact.close(); out_fungi.close(); out_other.close()


def get_query_clones_from_midwest_fasta(in_fasta):
    result = {} #{query: clone}
    for label, seq in MinimalFastaParser(in_fasta):
        query, desc = label.split(None, 1)
        clone = desc.split()[3]
        result[query] = clone
    return result


#########################
####run midwest
def main_midwest():
    midwest_fa = 'midwest.fa'
    midwest_blast = 'midwest.blast'
    midwest_bact_fa = 'midwest_bact.fa'
    midwest_fungi_fa = 'midwest_fungi.fa'
    midwest_other_fa = 'midwest_other.fa'

    """
    run_blast_all(['midwest.fa'], ['midwest_.blast'])
    set_trace()
    run_dump_nearest_common_lineages(file('midwest.fa'), file('midwest.blast'))
    """
    run_make_bact_fungi_other_fasta(file('midwest.fa'),
            file(midwest_bact_fa, 'w'), file(midwest_fungi_fa, 'w'),
            file(midwest_other_fa, 'w'))

midwest_clone_lineages_dump = 'midwest_clone_lineages.dump'

#########################
####run succession dataset
def main_succession():
    """
    Using: seqid_lineage_dump.
    """
    working_dir = '/home/zongzhi/Working/noah/succession'
    fasta = '%s/succession.fa' % working_dir
    blast = '%s/succession.blast' % working_dir
    lineage = '%s/succession.lineage' % working_dir

    bact_fa = '%s/succession_bact.fa' % working_dir
    fungi_fa = '%s/succession_fungi.fa' % working_dir
    other_fa = '%s/succession_other.fa' % working_dir
    bad_fa = '%s/succession_bad.fa' % working_dir

    seqid_lineages = dict(load(file(seqid_lineages_dump)))
    #blast and taxonomy
    #run_blast_all([fasta], [blast])
    nearest_common_lineages = nearest_common_lineages_from_blast(
            file(blast), seqid_lineages)

    #write lineages to a tab delimited file.
    out_lineage = file(lineage, 'w')
    out_lineage.write(
            '#seq_name\tcommon_lineage_of_the_nearest_blast_neighbors\n')
    for query, lineage in nearest_common_lineages:
        out_lineage.write('%s\t%s\n' % (query, '\t'.join(lineage)))
    out_lineage.close()

    nearest_common_lineages = dict(nearest_common_lineages)
    #write to seperate fasta files as requested.
    out_bact, out_fungi, out_other, out_bad = [file(fa, 'w')
        for fa in (bact_fa, fungi_fa, other_fa, bad_fa)]
    bad_queries = []
    for label, seq in MinimalFastaParser(file(fasta)):
        query, desc = label.split(None, 1)
        try:
            lin = nearest_common_lineages[query]
        except KeyError:
            no_hit_msg = '!!! no hit for  "-e 1e-10"'
            fasta = '>%s %s\n%s\n' % (query, no_hit_msg, seq)
            out_bad.write(fasta)
            continue

        fasta = '>%s %s\n%s\n' % (query, '; '.join(lin), seq)
        if lin[0:1] == ['Bacteria']:
            out_bact.write(fasta)
        elif lin[1:2] == ['Fungi']:
            out_fungi.write(fasta)
        else:
            print lin
            out_other.write(fasta)
    map(file.close, (out_bact, out_fungi, out_other, out_bad))
    print bad_queries
    set_trace()


######
# main
def main():
    #main_succession()
    return
    '''
    main_midwest()
    set_trace()

    seqid_lineages = dict(load(file(seqid_lineages_dump)))
    res = nearest_common_lineages_from_blast(
            file('%s/out_fasta/new_orlean.blast' % working_dir),
            seqid_lineages, threshold=0.9)
    res = [[name] + taxon for name, taxon in res]
    out_file = file('%s/query_new_orlean.xls' % working_dir, 'w')
    for row in res:
        out_file.write('%s\n' % '\t'.join(row))
    out_file.close()

    set_trace()

    t = load(file(query_lineages_dump))
    set_trace()
    #run_dump_seqid_lineages()
    #run_dump_query_toplineages()
    #run_dump_query_lineages()
    
    #query_lineages = load(file(query_lineages_dump))
    #run_query_common_taxon(query_lineages,
    #        'query_common_1.xls', threshold=1)
    #run_query_common_taxon(query_lineages,
    #        'query_common_0.9.xls', threshold=0.9)
    #run_query_common_taxon(query_lineages,
    #        'query_common_0.6.xls', threshold=0.6)

    query_toplineages = load(file(query_toplineages_dump))
    run_query_common_taxon(query_toplineages,
            'query_common_top_1.xls', threshold=1)
    run_query_common_taxon(query_toplineages,
            'query_common_top_0.9.xls', threshold=0.9)
    run_query_common_taxon(query_toplineages,
            'query_common_top_0.6.xls', threshold=0.6)
    set_trace()
    set_trace()
    '''
    fasta_path = working_dir + '/fasta/succession.fas'
    blast_path = working_dir + '/blast/succession.blast'
    set_trace()
    run_blast_all([fasta_path], [blast_path])
    set_trace()
    set_trace()
    t = TaxonomyFromBlast(
            BlastLines=file('blast/G910.out'),
            TaxidLineages=load('taxa.dump').items(),
            )
    #fasta_gis = list(t.eachQueryGi())
    #gi_taxids = list(t.eachGiTaxid())
    fasta_lineages = array(t.eachQueryLineage(), object)
    dump(t.GiTaxids, file('GiTaxids', 'w'))
    dump(t.QueryLineages, file('QueryLineages'))
    print fasta_lineages


if __name__ == '__main__':
    main()
    
