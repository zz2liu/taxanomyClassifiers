from get_rdp import *
from py_util.unit_test import TestCase, main, hook_pdb
from py_util.path_ import Path
hook_pdb()

class Tests(TestCase):
    def test_timetest(self):
        files = Path('data').glob('*.fa')
        out_files = [n+'.rdp' for n in files]
        res = timetest(files, out_files)
        self.p(res)

if __name__ == '__main__':
    main()
