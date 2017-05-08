"""rdp_classifier.py
"""
from cPickle import dump, load
from cogent.parse.fasta import MinimalFastaParser
from naive_classifier import NaiveClassifier


def rdp_fasta_parser(fasta):
    """yield each (seq as string, line age to genus as string).
        
    - fasta: fasta lines with rdp lineage to genus.
      ex: '>AB001783|S000020570 Root;Bacteria;Chlamydiae;Chlamydiae;'
      'Chlamydiales;Chlamydiaceae;Chlamydophila\n'
      'ATCCG\n'
    """
    for label, seq in MinimalFastaParser(fasta):
        try:
            name, lineage = label.split(None, 1)
        except ValueError, err:
            raise ValueError('label must have >= one whites, got %s' % label)
        yield seq, lineage


class RdpClassifier(NaiveClassifier):
    """train with fasta file and report result with lineage and bootstrap
    values.
    """
    def train(self, seqs, parser=rdp_fasta_parser, **train_kw):
        """setup ._word_posteriors.

        - seqs: input seq lines.
        - parser: a seq parser(seqs) -> each seq as str, genus as str
        **train_kw: pass to super.train.
        """
        super(RdpClassifier, self).train(parser(seqs), **train_kw)

    def classify(self, seqs, parser=MinimalFastaParser, **super_kw):
        """take seqs as input, return each (seqname, (genus, bootstrap))

        - seqs: query sequences.
        - parser: function(seqs) -> each (seqname, seq as str).
        **super_kw: pass to super.classify(,repeats=100)
        """
        super_fun = super(RdpClassifier, self).classify
        for name, seq in parser(seqs):
            yield name, super_fun(seq, **super_kw)
       

#### cmdline interface
def run_train(TYPE_SEQ, MODEL):
    """train using TYPE_SEQ to generate MODEL."""
    rdp = RdpClassifier()
    rdp.train(file(TYPE_SEQ))
    dump(rdp.__dict__, file(MODEL, 'w'))

def run_classify(QUERY_SEQ, MODEL):
    """stdout the classify results of QUERY_SEQ using MODEL."""
    rdp = RdpClassifier()
    rdp.__dict__ = load(file(MODEL))
    for name, (genus, bootstrap) in rdp.classify(file(QUERY_SEQ)):
        name, ori_lineage = name.split(None, 1)
        print ('%(name)s:\n'
                '%(ori_lineage)s\n'
                '%(genus)s, %(bootstrap)s%%' % locals())

def main():
    from sys import argv
    run_functions = {'train': run_train,
            'classify': run_classify
            }
    try:
        RUN, SEQ, MODEL = argv[1:]
        run = run_functions[RUN]
    except ValueError, err:
        raise ValueError('Usage: prog.py RUN_FUNCTION SEQ_FILE MODEL_FILE.\n'
            '  RUN_FUNCTION: one of %s.\n'
            '  SEQ_FILE: training sequences or query sequences.\n'
            '  MODEL_FILE: the model to build or using.'
            % (run_functions.keys(),))
    except KeyError, err:
        raise ValueError('RUN_FUNCTION must be one of %s, got %s'
                % (run_functions.keys(), RUN))

    run(SEQ, MODEL)

        
if __name__ == '__main__':
    main()




                

