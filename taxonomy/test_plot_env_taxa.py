"""test_plot_env_taxa.py

"""
from plot_env_taxa import *
from py_util.unit_test import TestCase, main, set_trace
from py_util.unit_test_plot import PlotTestCase
from numpy import arange, array
from textwrap import dedent
import pylab

class Test_get_name(TestCase):
    """test the support functions of get_name from lineage."""
    def setUp(self):
        self.certain = 'd p c o f'.split()

    def test_get_unclassified_name(self):
        self.assertEqual(get_unclassified_name(self.certain, 4),
                'f') # 4 family inx
        self.assertEqual(get_unclassified_name(self.certain[:-1], 4),
                'Unclassified_order_o')
        self.assertEqual(get_unclassified_name(self.certain[:-2], 4),
                'Unclassified_class_c')
        self.assertEqual(get_unclassified_name(self.certain[0:1], 4),
                'Unclassified_domain_d')
        self.assertEqual(get_unclassified_name([], 4),
                'Unclassified__')

        
        
class TestTitleEnvsI(TestCase):
    def setUp(self):
        self.tip_taxa = {
            'a': 'ABC AB A'.split(),
            'b': 'ABC AB B'.split(),
            'c': 'ABC CD C'.split(),
            'd': 'ABC CD'.split(),
            }
        self.env_lines = dedent("""\
            a\te1\t1
            b\te2\t2
            c\te3\t1
            d\te3\t2
            """).splitlines()
        self.phylum_env_counts = {
            'AB': {'e1': 1, 'e2': 2},
            'CD': {'e3': 3},
            }
        self.class_env_counts = {
            'A': {'e1': 1},
            'B': {'e2': 2},
            'C': {'e3': 1},
            '_CD': {'e3': 2},
            }


class Test_calc_title_env_matrix(TestTitleEnvsI):
    def setUp(self):
        super(Test_calc_title_env_matrix, self).setUp()
            
    def test_basic(self):
        tip_taxa, env_lines = self.tip_taxa, self.env_lines
        res = calc_title_env_matrix(tip_taxa, env_lines, 'phylum', )
        self.assert_equal(res, ([[1, 2, 0], [0, 0, 3]],
                ['AB', 'CD'], ['e1', 'e2', 'e3']))
        
        
class Test_calc_title_env_counts:#(TestTitleEnvsI, PlotTestCase):
    def setUp(self):
        super(Test_calc_title_env_counts, self).setUp()
            
    def test_basic(self):
        tip_taxa, env_lines = self.tip_taxa, self.env_lines
        res = calc_title_env_counts(tip_taxa, 'phylum', env_lines)
        self.assertEqual(res, self.phylum_env_counts)
        
        res = calc_title_env_counts(tip_taxa, 'class', env_lines)
        self.assertEqual(res, self.class_env_counts)

    ## break for now
    def test_fill_title_envs(self):
        title_env_counts = self.phylum_env_counts
        result = fill_title_envs(title_env_counts,
                title_order=['CD', 'AB'],
                env_order='e2 e1 e3'.split())
        self.outputFig()
        
        result = fill_title_envs(title_env_counts)
        self.outputFig()

class Test_dist_title_envs(TestTitleEnvsI):
    def setUp(self):
        super(Test_dist_title_envs, self).setUp()
        
    def test_dist_title_envs(self):
        result = dist_title_envs(self.phylum_env_counts,
                title_order=['CD', 'AB'], env_order='e2 e1 e3'.split())
        self.p(result)




class Test_reorder(TestCase):
    def setUp(self):
        self.mat =  array([[11, 12,  0],
                [21,  0, 1.3],
                [ 0, 3.2,  0]])
        self.row_names=['1', '2', '3']
        self.col_names = ['a', 'b', 'c']
     
    def test_reorder_cols(self):
        res, rnames = reorder_cols(self.mat, self.col_names, 'mean')
        self.assert_equal(res, self.mat[:, [2,1,0]])
        self.assert_equal(rnames, ['c', 'b', 'a'])

               
    def test_reorder_rows(self):
        res, rnames = reorder_rows(self.mat, self.row_names, 'mean')
        self.assert_equal(rnames, ['3', '2', '1'])

        #test another order method
        res, rnames = reorder_rows(self.mat, self.row_names, 'max')
        self.assert_equal(rnames, ['3', '1', '2'])

        #test no rownames
        res, rnames = reorder_rows(self.mat, None, 'max')
        self.assert_equal(res, self.mat[[2,0,1]])
        self.assert_equal(rnames, None)

        #test reverse
        res, rnames = reorder_rows(self.mat, self.row_names, 'max', reverse=True)
        self.assert_equal(res, self.mat[[1,0,2]])
        self.assert_equal(rnames, ['3', '1', '2'][::-1])
        
class Test_merge_rows(TestCase):
    def setUp(self):
        self.mat = array([[11, 12,  0],
                [21,  0, 23],
                [ 0, 32,  0]])
        self.row_names=['1', '2', '3']

    def test_basic(self):
        a2d, rownames = self.mat, self.row_names
        # first two merged
        mat, rnames = merge_rows(a2d, rownames=rownames, merge_threshold=29)
        self.assert_equal(mat,  [[ 0, 32,  0], [32, 12, 23]])
        self.assert_equal(rnames, ['3','1+2'])
                
    def test_first_one_merged(self):
        a2d, rownames = self.mat, self.row_names
        # first one merged, all elem of row < thres
        res = merge_rows(a2d, rownames=rownames, merge_threshold=20,
                merged_name='merged')
        self.assert_equal(res, (a2d[[1,2,0]], ['2','3','merged']))
    
    def test_no_merged(self):
        a2d, rownames = self.mat, self.row_names
        # no row merged
        res = merge_rows(a2d, rownames=rownames, merge_threshold=10)
        self.assert_equal(res, (a2d, rownames))

    def test_all_merged(self):
        a2d, rownames = self.mat, self.row_names        
        # all rows merged
        res = merge_rows(a2d, rownames=rownames, merge_threshold=39)
        self.assert_equal(res, ([a2d.sum(0)], ['1+2+3']))


    def test_13_merged(self):
        a2d = array([[11, 12,  0],
                [21,  0, 23],
                [ 0, .32,  0]])
        rownames=['1', '2', '3']
        # 1,3 merged
        res = merge_rows(a2d, rownames=rownames, merge_threshold=19)
        self.assert_equal(res, (
                [[ 21.  ,   0.  ,  23.  ],
                 [ 11.  ,  12.32,   0.  ]],
                [ '2', '1+3'])
                )

class Test_cum_cols(TestCase):
    def test_basic(self):
        m = arange(9).reshape(3,3).astype(float)
        exp = [[  0.,   0.,   0.,],
                [  0.,   1.,   2.,],
                [  3.,  5.,   7.,],
                [  9.,  12.,  15.,],]
        self.assert_equal(cum_cols(m), exp)
        self.assert_equal(cum_cols(m, add_base=False), exp[1:])

class Test_fill_plot:#(PlotTestCase): ## test a wrapper
    Debug = False
    def test_basic(self):
        a2d = array([[1.1, 1.2,  0],
                [2.1,  0, 2.3],
                [ 0, .32,  0]])
        rownames=['1', '2', '3']
        colnames=['a', 'b', 'c']
        fill_plot(a2d, rownames, colnames)
        self.outputFig()

        fill_plot(a2d, rownames, colnames, norm=False)
        self.outputFig()

        fill_plot(a2d, rownames, colnames, merge_threshold=1.0)
        self.outputFig()


class Test_fillplot_heatplot(PlotTestCase): ## test a wrapper
    Debug = False
    def setUp(self):
        self.a2d = array([[1.1, 1.2,  0],
                [2.1,  0, 2.3],
                [ 0, .32,  0]])
        self.rownames=['1', '2', '3']
        self.colnames=['a', 'b', 'c']
    
    def test_basic(self):
        fillplot(self.a2d, self.rownames, self.colnames)
        self.outputFig()
        fillplot(self.a2d)
        self.outputFig()
    
    def test_heatplot(self):
        heatplot(self.a2d, self.rownames, self.colnames)
        self.outputFig()
        heatplot(self.a2d)
        self.outputFig()

    def test_fillplot_matrix(self):
        #self.Debug=True
        fillplot_matrix(self.a2d, self.rownames, self.colnames)
        self.outputFig()
        fillplot_matrix(self.a2d)
        self.outputFig()

if __name__ == '__main__':
    main()
