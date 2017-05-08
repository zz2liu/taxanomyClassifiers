"""plot_env_taxa.py

09/26/07 created from 
todo: get rid of dep on RdpTaxonomy, Dict2D

todo: finish the mat + row_names + col_names version
todo: add method to TableLite - merge, reorder
"""
from __future__ import division
from functools import partial
from warnings import warn
from numpy import  array, asarray, zeros, concatenate, cumsum, iterable, argsort
import pylab
from cogent.util.dict2d import Dict2D

from taxonomy.util import TITLES
from unifrac.env_analysis_table import EnvAnalysis
from my_util.draw.plot import fill_between, DISTINCTIVE_COLORS
from my_util.maths.distances import distance_matrix_from_points

# tobe mvd to a more general folder
def matrix_from_dict2d(dict_2d, row_order_f=sorted, col_order_f=sorted):
    """return a 2darray from a dict2d.

    - dict2d: a dict of {rowname: {colname: value}}.
    - row_order_f=sorted: a fun(all_rownames)-> rownames in order
    - col_order_f=sorted: a fun(all_colnames) -> colnames in order
    """
    def get_colnames(dict_2d):
        result = set()
        for colname_value in dict_2d.values():
            result.update(colname_value.keys())
        return result

    rownames = dict_2d.keys()
    colnames = get_colnames(dict_2d)

    row_order = row_order_f(rownames)
    col_order = col_order_f(colnames)

    result = Dict2D(dict_2d,
            RowOrder=row_order, ColOrder=col_order,
            Default=0, Pad=True).toLists()
    return array(result), row_order, col_order


###
# support functions
def _lineage_lacks(taxon, idx):
    lineage = taxon[:idx+1]
    num_uncertain = idx+1 - len(lineage)
    return lineage, num_uncertain

def get_lineage_name(taxon, idx):
    lineage, num_uncertain = _lineage_lacks(taxon, idx)
    return '_'*num_uncertain + ';'.join(lineage)

def get_last_name(taxon, idx):
    lineage, num_uncertain = _lineage_lacks(taxon, idx)
    return '_'*num_uncertain + (lineage[-1] if lineage else '')

def get_title_name(taxon, idx):
    name = ''.join(taxon[idx:idx+1])
    if name in ('', 'incertae'):
        name = 'unknown_%s' % TITLES[idx]
    return name

def get_tip_titles(tip_taxa, title, get_name='last'):
    """return {tipname: title_name}.
    
    - tip_taxa: {tipname: taxon as a list of taxs}
    - title: title name or idx
    - get_name='last': 'last', 'lineage' or a fun(taxon, title_idx) -> title_name.
    """
    title_idx = (title if isinstance(title, int)
            else TITLES.index(title))
    get_name = (get_last_name  if get_name=='last'
            else get_lineage_name  if get_name=='lineage'
            else get_name)
    #{tip: title}
    tip_titles = {}
    for tip, taxon in tip_taxa.iteritems():
        tip_titles[tip] = get_name(taxon, title_idx)
    return tip_titles

## 26/07/07 mv here from run_taxonomy_from_arb.py
# tobe improved by take tip_titles and use default dict.
def calc_title_env_counts(tip_taxa, title, env_file,
        return_more=False, rescue_title=None, **kw):
    """return a dict of {title: {env: count}}

    - tip_taxa: {tip: taxon as str list}
    - title: a str like 'class', 'family', and so on.
    - env_file: envfile used for unifrac, warn if tipname not in env_file.
    - rescue_title: an optional function to rescue some unknown titles.
    Deps: RdpTaxonomy, EnvAnalysis
    return_more: for test purpose only.
    """
    
    tip_titles = get_tip_titles(tip_taxa, title, **kw)

    #{tip: [(env, count)]}
    tip_env_counts = EnvAnalysis(envs=env_file).seqEnvs()

    #{title: [(env, count),]} -- with possible env duplicates
    title_env_counts = {}
    bad_tips = []
    for tip, title in tip_titles.items():
        env_counts = tip_env_counts.get(tip, None)
        if not env_counts:
            bad_tips.append(tip)
            continue
        title_env_counts.setdefault(title, []).extend(env_counts)
    if bad_tips:
        warn('tipnames %s not found in env_file' % bad_tips)
            
    # {title: {env: count}} -- with counts merged
    result = {} 
    for title, env_counts in title_env_counts.items():
        val = {} # {env: count,}
        for env, count in env_counts:
            val[env] = val.get(env, 0) + count
        result[title] = val

    if return_more:
        return result, {'tip_titles': tip_titles,
                'tip_env_counts': tip_env_counts,
                'title_env_counts': title_env_counts}
    return result





def norm_cols(data):
    """return a float 2darray with each col normalized by /sum that col_order_f
    """
    return data / data.sum(0)
#norm_methods = {'to01': _norm_to01}
 
def merge_rows(data, rownames=None, merge_threshold=0.05, merged_name=None):
    """return a fillplot after norm and merging.

    - data: a float 2darray, no check.
    - rownames: rownames in order for the matrix, no check.
    - merge_threshold: a row will be merged if all its values are below the
     threshold.
    - merged_name: the row_name of the merged row, default to be
     '+'.join(row_names merged)

    Note: the merged_row is put at the bottom.
    """
    rownames_ori = rownames
    if rownames is None:
        rownames = ['']*len(data)
 
    major_rows = []; minor_rows = []
    major_names = []; minor_names = []
    for row, rowname in zip(data, rownames):
        if all(row < merge_threshold):
            minor_rows.append(row)
            minor_names.append(rowname)
        else:
            major_rows.append(row)
            major_names.append(rowname)
    if not minor_rows:
        return data, rownames
    
    merged_row = array(minor_rows).sum(0).tolist()
    merged_name = merged_name or '+'.join(minor_names)
    result_rows =  major_rows + [merged_row]
    rownames = major_names + [merged_name]
    if rownames_ori is None: rownames = None
    return array(result_rows), rownames

def cum_cols(data, add_base=True):
    """return a data.cumsum(0) with optional baseline."""
    cum_data = data.cumsum(axis=0)
    if add_base:
        baseline = zeros(cum_data.shape[1])
        cum_data = concatenate([[baseline], cum_data], axis=0)
    return cum_data

def reorder_rows(mat, row_names=None, orderby='mean', reverse=False):
    """reorder the rows of a mat return a new mat with names.

    - mat, row_names: a matrix with its row_names
    - orderby='mean': one of [mean, std, min, max] or a seq of names.
    """
    if isinstance(orderby, str):
        assert orderby in ['mean', 'std', 'min', 'max', 'sum', 'var']
        method = getattr(mat, orderby)
        rule = method(1)
        idx_order = argsort(rule)
    else: # orderby is a list of names
        try:
            idx_order = [row_names.index(n) for n in orderby]
        except (IndexError, AttributeError), e:
            raise ValueError('name %r not found in names %s' %
                    (n, names))

    if reverse: idx_order = idx_order[::-1]
    mat = mat.take(idx_order, axis=0)
    if row_names:
        row_names = asarray(row_names)[idx_order]
    return mat, row_names

def reorder_cols(mat, col_names=None, orderby='mean', reverse=False):
    """reorder the cols of a mat return a new mat with names.

    - mat, col_names: a matrix with its col_names
    - orderby='mean': one of [mean, std, min, max] or a seq of names.
    """
    if isinstance(orderby, str):
        assert orderby in ['mean', 'std', 'min', 'max', 'sum', 'var']
        method = getattr(mat, orderby)
        rule = method(0)
        idx_order = argsort(rule)
    else: # orderby is a list of names
        try:
            idx_order = [col_names.index(n) for n in orderby]
        except (IndexError, AttributeError), e:
            raise ValueError('name %r not found in names %s' %
                    (n, col_names))

    if reverse: idx_order = idx_order[::-1]
    mat = mat.take(idx_order, axis=1)
    if col_names:
        col_names = asarray(col_names)[idx_order]
    return mat, col_names


def fillplot(mat, row_names=None, col_names=None, **kw):
    """fillplot each colum of a matrix."""
    cum_mat = cum_cols(mat, add_base=True)
    return fill_between(cum_mat, row_names=row_names, col_names=col_names,
            **kw)

def heatplot(mat, row_names=None, col_names=None, **kw):
    """heatplot a matrix."""
    result = [pylab.matshow(mat, **kw)]

    ## set x, y ticks by col, row names.
    if row_names:
        result.append(pylab.yticks(range(mat.shape[0]), row_names))
    if col_names:
        result.append(pylab.xticks(range(mat.shape[1]), col_names))
    return result


####
# the out most wrappers

def calc_title_env_matrix(tip_taxa, env_file, level='class', get_name='last'):
    """return a (title by env matrix of counts, row_names, col_names).

    - tip_taxa: {tipname: [lineage]}
    - env_file: a unifrac env file, each line as 'seqname\tenvname\tcount\n'
    - level: phylum, class, order, family or genus
    - get_name='last': a method to get a name from a lineage, tobe rm.
    """
    title_env_counts = calc_title_env_counts(tip_taxa, level, env_file, get_name=get_name)
    return matrix_from_dict2d(title_env_counts)

            
def fillplot_matrix(mat, row_names=None, col_names=None,
        row_order=('mean', True), col_order=None,
        norm=True, merge=(0.02, 'Others'),
        **fillplot_kw):
    """fillplot a title by env matrix of counts.
    Return (mat, row_names, col_names) for later use.

    - mat, row_names, col_names: the matrix with names
    - row_order=('mean', True): pass to reorder_rows(,orderby, reverse)
    - col_order=None: reorder method of cols.
    - norm=True: normalize each col to 0..1
    - merge=(0.02, 'Others'): the merge threshold and merged_name for rows.
    - **fillplot_kw: pass to fillplot(,colors=, xticks_pars=, legend_pars=,...)
    
    Note: The order of mat manipulations is: reorder->norm->merge
    """
    if row_order:
        mat, row_names = reorder_rows(mat, row_names, *row_order)
    if col_order:
        mat, col_names = reorder_cols(mat, col_names, *col_order)
    if norm:
        mat = norm_cols(mat)
    if merge:
        mat, row_names = merge_rows(mat, row_names, *merge)

    fillplot(mat, row_names, col_names, **fillplot_kw)
    return mat, row_names, col_names
    
        
    














    
### deprecated functions below
def fill_title_envs(title_env_counts, title_order=None, env_order=None,
        colors=DISTINCTIVE_COLORS,
        title_labels=None, env_labels=None,
        norm=True, verbose=False, **fill_plot_kw):
    """fill the polygons between each env_counts accumulated.

    - title_env_counts: {title: {env: count}}
    - title_order, env_order: order of the title_keys, env_keys.
    - title_labels=None, env_labels=None: provide if you need labels different
      from keys.
    - **fill_plot_kw: (norm=, addbase=, xticks_pars=, lengend_pars=,  )
    """ 
    data, title_order, env_order = get_title_env_matrix(
            title_env_counts, title_order=title_order, env_order=env_order)
    title_labels = title_labels or title_order
    env_labels = env_labels or env_order

    #fill polygons
    fill_plot(data, rownames=title_labels,
            colnames=env_labels, colors=colors, **fill_plot_kw)


def dist_title_envs(title_env_counts, title_order, env_order):
    """return a distance matrix from {title: {env: count}}.

    Percent abundance of each title as env point coordinates.
    - title_order, env_order: a seq in order.
    """
    if not iterable(title_order): raise TypeError()
    if not iterable(env_order): raise TypeError()
    
    env_cols = get_title_env_matrix(title_env_counts,
            title_order=title_order, env_order=env_order)[0]
    env_cols = array(env_cols, float)
    
    env_cols =  env_cols/env_cols.sum(0)
    #points = title_env_count_matrix.T #rows as points
    dist_m = distance_matrix_from_points(env_cols, bycol=True)
    return dist_m

    
            
def fill_plot(data, rownames, colnames, norm=True,
        merge_threshold=None, merged_name=None, **fill_between_kw):
    """a wrapper to fill_plot.
    
    - data, rownames, colnames: a float table.
    - norm=True:
    - merged_threshold=None:
    """
    if norm:
        data = norm_cols(data)
    if merge_threshold is not None:
        data, rownames = merge_rows(data, rownames,
                merge_threshold=merge_threshold,
                merged_name=merged_name)
    cum_data = cum_cols(data, add_base=True)
    return fill_between(cum_data, row_names=rownames, col_names=colnames,
            **fill_between_kw)
                
####
# special funs        
        
sorted_reverse = partial(sorted, reverse=True)
def get_title_env_matrix(title_env_counts, title_order=sorted_reverse,
        env_order=sorted):
    """return (title_by_env 2darray, title_orer, env_order).

    - title_env_counts: a dict of {title: {env: count}}.
    - title_order=sorted_reverse: a fun(titles) -> title_order
      or a seq of title_order.
    - env_order=sorted: a fun(envnames)->env_order or a seq of env_order
    """
    row_order_f = (title_order if callable(title_order)
            else lambda x: title_order)
    col_order_f = (env_order if callable(env_order)
            else lambda x: env_order)
    return matrix_from_dict2d(title_env_counts,
            row_order_f=row_order_f, col_order_f=col_order_f)    
                
def run(TIPTAXA, ):
    pass

    
    
