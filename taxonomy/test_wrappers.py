"""test_wrappers.py
"""
from wrappers import *
from py_util.unit_test import TestCase, main, set_trace

class Test_dump_tip_taxa_from_groupnames(TestCase):
    def test_basic(self):
        phylo_tree = '/home/zongzhi/Projects/misc/zongzhi/taxonomy/test_data/test_from_groupnames.ntree'
        tip_taxa_pickle = './~tmp'
        result = dump_tip_taxa_from_groupnames(phylo_tree, tip_taxa_pickle)
        self.assertEqual(len(result), 24)
        result = dump_tip_taxa_from_groupnames(phylo_tree, tip_taxa_pickle,
                selector=lambda n: n.Name.startswith('VA1'))
        self.assertEqual(len(result), 15)


if __name__ == '__main__':
    main()
