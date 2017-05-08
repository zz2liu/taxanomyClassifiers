from py_util.unit_test import TestCase, main, set_trace
from gi_to_taxid import TaxidFinder, create_db

from functools import partial
import os

test_db_path = 'gi_taxid_test.db'
try: os.remove(test_db_path)
except OSError: pass
create_db(test_db_path, 'gi_taxid_test.dmp')
class Test_TaxidFinder(TestCase):
    def test_basic(self):
        tf = TaxidFinder(test_db_path)
        self.assertEqual(tf(gis_test), result_test)
        self.assertEqual(tf(gis_test+[0]), result_test) #3 is ignored
        self.assertEqual(tf([]), [])

    def test_bad_gis(self):
        tf = TaxidFinder(test_db_path)
        self.assertEqual(tf(gis_test, return_more=True), (result_test, []))
        self.assertEqual(tf(gis_with_bad_test, return_more=True), 
                (result_test, [0]))
        self.assertEqual(tf([]), [])



### test data
gis_test = [2, 323587, 708579]
gis_with_bad_test = [0, 2, 323587, 708579] #o not in table
result_test = [(2, 9913), (323587, 11053), (708579, 9606)]
if __name__ == '__main__':
    #print gis_to_taxids(gis_test, 'gi_taxid.db')
    main()
