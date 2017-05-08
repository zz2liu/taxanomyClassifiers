from py_util.unit_test import TestCase, main, set_trace
from taxid_to_lineage import (sqlite3, create_db, LineageFinder)

from py_util import TempFileName
from _test_taxid_to_lineage import nodes_test, names_test

test_db_name = TempFileName(prefix='TempFileName_test_taxid_to_lineage_')
db_path = test_db_name
create_db(nodes_test, names_test, test_db_name)

def print_first_rows(db_path, table_name, n=5):
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute("select * from %s" % table_name)
    for i, row in enumerate(cur):
        if i > n: break
        print row
    con.close()

class test_create_db(TestCase):
    def test_create_db(self):
        print_first_rows(test_db_name, 'nodes', 20)

class Test_LineageFinder(TestCase):
    def setUp(self):
        self.finder = LineageFinder(db_path)
        self.cur = self.finder._cur

    def test___init(self):
        self.p(self.finder._cur, self.finder._con)
        
    def test__get_parent_id(self):
        cur = self.cur
        self.p(LineageFinder._get_parent_id(cur, 18))

    def test_get_ancestors(self):
        self.p(LineageFinder._get_ancestor_ids(self.cur, 18))

    def test_get_id_rank_names(self):
        cur = self.cur
        result = LineageFinder._get_id_rank_names(cur, [18, 10, 2, 1])
        self.p(result)

    def test_taxids_to_lineages(self):
        res = self.finder.findLineages([18, 10, 2])
        self.p(list(res))

    def test_taxid_to_lineage(self):
        res = LineageFinder._taxid_to_lineage(self.cur, 18)
        self.p(res)

if __name__ == '__main__':
    #LineageFinder.createDb('taxdump/nodes.dmp', 'taxdump/names.dmp',
    #        'taxdump/node_names.db')

    #finder = LineageFinder('taxdump/node_names.db')
    #res = finder([451511, 307270, 235235, 115162, 12122, 6, 2, 1])
    #for e in res:
    #    print e
    
    main()
