from my_util.unit_test import TestCase, main, set_trace
from taxid_to_lineage import (sqlite3, create_db, taxid_to_lineage,
        taxids_to_lineages, _get_parent_id, get_ancestor_ids, get_id_rank_names)

from my_util import TempFileName
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

class Test_taxids_to_lineages(TestCase):
    def setUp(self):
        self.con = sqlite3.connect(db_path)
        self.cur = self.con.cursor()

    def tearUp(self):
        self.con.close()

    def test_create_db(self):
        print_first_rows(test_db_name, 'nodes', 20)
        
    def test__get_parent_id(self):
        cur = self.cur
        self.p(_get_parent_id(cur, 18))

    def test_get_ancestors(self):
        self.p(get_ancestor_ids(self.cur, 18)) #, test_db_name))

    def test_get_id_rank_names(self):
        cur = self.cur
        result = get_id_rank_names(cur, [18, 10, 2, 1])
        self.p(result)

    def test_taxids_to_lineages(self):
        res = taxids_to_lineages([18, 10, 2], test_db_name)
        self.p(list(res))

    def test_taxid_to_lineage(self):
        res = taxid_to_lineage(self.cur, 18)
        self.p(res)

if __name__ == '__main__':
    main()
