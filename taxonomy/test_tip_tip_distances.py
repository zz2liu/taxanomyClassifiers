from tip_tip_distances import (get_tip_tip_distances_dict,
        get_tip_tip_distances_array, get_tip_tip_distances)
from py_util.unit_test import TestCase, main
from StringIO import StringIO
from py_util.util import set_trace
from py_util.tree_new import DndParser

class Test_tip_tip_distances_I(object):
    def setUp(self):
        self.root_std = DndParser(StringIO(tree_std))
        self.root_one_level = DndParser(StringIO(tree_one_level))
        self.root_two_level = DndParser(StringIO(tree_two_level))
        self.root_one_child = DndParser(StringIO(tree_one_child))
        #set_trace()

    def test_one_level(self):
        res = self.fun(self.root_one_level)
        self.pp(res)

    def test_two_level(self):
        res = self.fun(self.root_two_level)
        self.pp(res)

    def test_one_child(self):
        res = self.fun(self.root_one_child)
        self.pp(res)

    def test_std(self):
        res = self.fun(self.root_std)
        self.pp(res)


class Test_tip_tip_distances_dict(Test_tip_tip_distances_I, TestCase):
    def setUp(self):
        self.fun = get_tip_tip_distances_dict
        super(Test_tip_tip_distances_dict, self).setUp()

class Test_tip_tip_distances_array(Test_tip_tip_distances_I, TestCase):
    def setUp(self):
        self.fun = get_tip_tip_distances_array
        super(Test_tip_tip_distances_array, self).setUp()

    def test_std(self):
        dist, tips = self.fun(self.root_std)
        tips = [tip.Name for tip in tips]
        self.assertEqual(dist, tree_std_dist)
        self.assertEqual(tips, tree_std_tips)

    def test_one_child(self):
        dist, tips = self.fun(self.root_one_child)
        tips = [tip.Name for tip in tips]
        self.assertEqual(dist, tree_one_child_dist)
        self.assertEqual(tips, tree_one_child_tips)

class Test_tip_tip_distances(Test_tip_tip_distances_I, TestCase):
    def setUp(self):
        self.fun = get_tip_tip_distances
        super(Test_tip_tip_distances, self).setUp()

    def test_std(self):
        dist, tips = self.fun(self.root_std, verbose=True)
        tips = [tip.Name for tip in tips]
        self.assertEqual(dist, tree_std_dist)
        self.assertEqual(tips, tree_std_tips)

    def test_one_child(self):
        dist, tips = self.fun(self.root_one_child)
        tips = [tip.Name for tip in tips]
        self.assertEqual(dist, tree_one_child_dist)
        self.assertEqual(tips, tree_one_child_tips)

### test data
tree_std = """\
        ((a:1, b:2, c:3)abc:0.1, (d:4, (e:5, f:6)ef:0.2)def:0.3);
"""
tree_std_dist = \
      [[  0. ,   3. ,   4. ,   5.4,   6.6,   7.6],
       [  3. ,   0. ,   5. ,   6.4,   7.6,   8.6],
       [  4. ,   5. ,   0. ,   7.4,   8.6,   9.6],
       [  5.4,   6.4,   7.4,   0. ,   9.2,  10.2],
       [  6.6,   7.6,   8.6,   9.2,   0. ,  11. ],
       [  7.6,   8.6,   9.6,  10.2,  11. ,   0. ]]
tree_std_tips = ['a', 'b', 'c', 'd', 'e', 'f']

tree_one_level = """(a:1, b:2, c:3)abc;"""

tree_two_level = """((a:1, b:2, c:3)abc:0.1, d:0.3)abcd;"""

tree_one_child = """((a:1, b:2, c:3)abc:0.1, (d:0.2)d_:0.3)abcd;"""
tree_one_child_dist = \
      [[ 0. ,  3. ,  4. ,  1.6],
       [ 3. ,  0. ,  5. ,  2.6],
       [ 4. ,  5. ,  0. ,  3.6],
       [ 1.6,  2.6,  3.6,  0. ]]
tree_one_child_tips = ['a', 'b', 'c', 'd']


if __name__ == '__main__':
    main()
