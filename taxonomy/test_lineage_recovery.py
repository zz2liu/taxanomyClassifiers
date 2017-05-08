from __future__ import division
from lineage_recovery import *
from py_util.unit_test import TestCase, main, set_trace
from numpy import *

class MiscTests(TestCase):
    def test_get_lineage_recovery(self):
        inp_exp = [ #(obs_lin, exp_lin), exp
                ( (['a', 'bb', 'c'], ['a', 'b', 'c']),  [1, 0, 1] ),
                ( (['a', 'bb', None], ['a', 'b', 'c']), [1, 0, -1] ),
                ( (['a', 'bb', 'c'], ['a', None, 'c']), [1, -1, 1] ),
                ]
        for inp, exp in inp_exp:
            self.assertEqual(array(get_lineage_recovery(*inp), int), exp)

    def test_calc_recovery_scores(self):
        res = calc_recovery_scores( [[1, 0, 1],
                                     [1, 0, -1],
                                     [1, -1, 0]])
        self.assertEqual(res,
                [(1.0, 0.0, 0.5), (1, 2/3, 2/3)])


if __name__ == '__main__':

    main()
