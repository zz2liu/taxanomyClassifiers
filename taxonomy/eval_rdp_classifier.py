"""eval_rdp_classifier.py
"""
def parse_names_and_seq_genera(seqs):
    """return [seq_name] and [(seq, genus),] in the same order.
    """
    names, seq_lineages = [], []
    for label, seq in MinimalFastaParser(seqs):
        name, lineage = label.split(None, 1)
        names.append(name)
        seq_lineages.append((seq, lineage))
    return names, seq_lineages

class RdpLeaveOneOutEvaluation(object):
    """evaluation class.
    """
    def __init__(self, seqs):
        """seqs must be a list of seqs"""
        self._names, self._seq_genera = parse_names_and_seq_genera(seqs)
        self.classifier = NaiveClassifier()

    def eval(self):
        names, seq_genera, classifier = \
                self._names, self._seq_genera, self.classifier
        for i in range(len(names)):
            classifer.train(seq_genera[:i] + seq_genera[i+1:])
            name = names[i]; seq, genus = seq_genera[i]
            genus_, bootstrap = classifer.classify(seq)
            print ('%(name)s: %(genus)s\n'
                    '%(bootstrap)s - %(genus_)s\n')


