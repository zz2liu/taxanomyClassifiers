"""create dpcofg NameArray tab file from file of [(seqname, dpcofg as list),]

Note: run_dpcofg_list_to_plots.py which has an upper level wrapper.


Use create_dpcofg_table_files(..., truncate=False) to count the names lost
after clipping.

each clip_name is expanded to >=1 names to match those before clipping.
each lineage is fixed by padding the missing tax in the middle with its
ancestor name + '_'.

03/04/08 mv here from working, disable some highly coupled functions

todo: get rid of hard code, add test code.
"""
from warnings import warn
from pdb import set_trace
from pickle import load, dump
from numpy import disp, zeros
from pdb import set_trace
from os import path

from py_util.numpy_ import NamedArray
from taxonomy.lineage_tree import pad_mid_missing_tax
from util import TITLES as RANK_ORDER
#from util import ( RANK_ORDER, ORI_CLIP, CLIP_DICT_DIR,
        #method_dataset_clip_from_fname)
from py_util.unit_test import hook_pdb
hook_pdb()


def _clip_from_filename(filename):
    """return clip from a groupname-dpcofg_list filename."""
    return path.basename(filename).split('.')[0]

#another main horse
def create_dpcofg_table_files_using_clip_dict(ori_filename, clip_filenames,
        clip_uniques, out_suffix='.expanded_tab',
        clip_from_filename=_clip_from_filename, expand_ori=True, **kw):
    """create each dpcofg table file with the same order as ori_file.
    but it can have less rows.

    - ori_filename: the pickle file with ori classifications as
      [(seqname, dpcofg as a list)]
    - clip_filenames: the pickle files with clip classifications.
    - clip_uniques: a dict of {clip: {seqname: expanded names}}
    - out_suffix='.expanded_tab': used to create the output table filename.

    Using: create_ori_dpcofg_table, create_clip_dpcofg_table
    """
    kw.setdefault('verbose', True)
    kw.setdefault('strict', True)

    if kw['verbose']:
        disp(ori_filename + ' ...')
    ori_list = load(open(ori_filename))
    clip_dict = ( None if not expand_ori
            else clip_uniques[clip_from_filename(ori_filename)])
    ori_table = create_ori_dpcofg_table(ori_list, clip_dict)
    ori_table.toTabFile(ori_filename + out_suffix,
            header=True, lefter=True)

    #ori_table is now safe to be cell-modified.
    for clip_filename in clip_filenames:
        clip = clip_from_filename(clip_filename)
        if kw['verbose']:
            disp('%s, %s ...' % (clip_filename, clip))
        clip_list = load(open(clip_filename))
        try: clip_dict = clip_uniques[clip]
        except: set_trace()
        clip_table = create_clip_dpcofg_table(ori_table,
                clip_list, clip_dict, **kw)
        clip_table.toTabFile(clip_filename + out_suffix,
                header=True, lefter=True)

#main horse
def create_dpcofg_table_files(ori_filename, clip_filenames,
        clip_dict_filenames, ori_dict_filename=None,
        out_suffix='.expanded_tab', **kw):
    """create each dpcofg table file with the same order as ori_file.
    but it can have less rows.

    - ori_filename: the pickle file with ori classifications as
      [(seqname, dpcofg as a list)]
    - clip_filenames: the pickle files with clip classifications.
    - clip_dict_filenames: the pickle files with {seqname: seqnames represented
      as a list}
    - ori_dict_filename: if provided expand the ori_list first.
    - out_suffix='.expanded_tab': used to create the output table filename.
    **kw: pass to create_clip_dpcofg_table(,clip_dict=)

    Using: create_ori_dpcofg_table, create_clip_dpcofg_table
    """
    kw.setdefault('verbose', True)
    kw.setdefault('strict', True)

    if kw['verbose']:
        disp(ori_filename + ' ...')
    ori_list = load(open(ori_filename))
    ori_dict = ( load(open(ori_dict_filename)) if ori_dict_filename
            else {} )
    ori_table = create_ori_dpcofg_table(ori_list, ori_dict)
    ori_table.toTabFile(ori_filename + out_suffix,
            header=True, lefter=True)

    #ori_table is now safe to be cell-modified.
    for clip_filename, clip_dict_filename in zip(
            clip_filenames, clip_dict_filenames):
        if kw['verbose']:
            disp(clip_filename + ' ...')
        clip_list = load(open(clip_filename))
        clip_dict = load(open(clip_dict_filename))
        clip_table = create_clip_dpcofg_table(ori_table,
                clip_list, clip_dict, **kw)
        clip_table.toTabFile(clip_filename + out_suffix,
                header=True, lefter=True)

def create_ori_dpcofg_table(ori_list, clip_dict={}, field_names=RANK_ORDER):
    """return a dpcofg_table as 2dNameArray for original classifications.
    """
    # make the ori_table
    rec_names, data = [], []
    for name, dpcofg in ori_list:
        dpcofg = pad_mid_missing_tax(dpcofg)
        name_expanded = clip_dict.get(name, [name])
        rec_names.extend(name_expanded)
        data.extend([dpcofg] * len(name_expanded))
    ori_table = NamedArray(data, object, names=[rec_names, field_names])
    return ori_table

def create_ori_dpcofg_table_old(ori_list, field_names=RANK_ORDER):
    """return a dpcofg_table as 2dNameArray for original classifications.
    """
    # make the ori_table
    rec_names, data = [], []
    for name, dpcofg in ori_list:
        dpcofg = pad_mid_missing_tax(dpcofg)
        rec_names.append(name)
        data.append(dpcofg)
    ori_table = NamedArray(data, object, names=[rec_names, field_names])
    return ori_table

def create_clip_dpcofg_table(ori_table, clip_list, clip_dict,
        verbose=True, strict=True):
    """return a dpcofg_table (NamedArray) for clip classifications.
    the result table is with the same shape and names as that of ori_table.

    - strict=False: if true, do not truncate the clip_table.
    Note: ori_table cells will be changed in-place.
    """
    # actually all I need is the rownames of ori_table
    (nrow, ncol), names = ori_table.shape, ori_table._names
    result = [['']*ncol]*nrow
    result = NamedArray(result, object,  names=ori_table._names)

    for name, dpcofg in clip_list:
        dpcofg = pad_mid_missing_tax(dpcofg)
        try: names_expand = clip_dict[name]
        except KeyError, e: #many seqname becomes seqname! after clearcut
                #'SBYC466_clip_83_1326' -> 'SBYC466_clip_83_1326\x19'
            name = name.strip('!\x19')
            names_expand = clip_dict[name]
        try: result[names_expand] = dpcofg
        except KeyError, e: #some seqs not assigned as ori
            names_expand = list(set(names_expand) & set(names[0]))
            result[names_expand] = dpcofg
    if strict:
        return result
    else:
        warn('alway use strict.', DeprecationWarning)
        return result[all_names]

#def get_files_new(dataset, files, clip_dict_dir, fname_parser):
    #"""quick code.
    #eg: get_files('ruth', glob('*.dpcofg_list'))
    #"""
    #ori_file = None; clip_files=[]; dict_files = []
    #for name in files:
        #method, ds, clip = fname_parser(name)
        #if ds != dataset:
            #continue
        #if clip == ORI_CLIP:
            #ori_file = name
        #else:
            #clip_files.append(name)
            #dict_files.append('%s/%s__%s__unique.dict'
                    #% (clip_dict_dir, dataset, clip))
    #return ori_file, clip_files, dict_files

#def to_tables(files, fname_parser=method_dataset_clip_from_fname,
        #clip_dict_dir=CLIP_DICT_DIR, **kw):
    #"""
    #*a: parse to get_files(dataset, suffix, ori_between, clip_between_glob)
    #**kw: strict=True: count the names lost after clipping as a coverage loss.
    #"""
    #for dataset in ['ruth', 'relman', 'gn']:
        #ori_file, clip_files, dict_files = get_files_new(dataset, files,
                #fname_parser, clip_dict_dir)
                ##'dpcofg_list', '.fa.blast_.', '__*.')
        #create_dpcofg_table_files(ori_file, clip_files, dict_files,
                #**kw) #strict=strict)


#def main(argv):
    #"""usage eg:
    #call('python PROG.py to_tables \"%s\"' % glob('*.dpcofg_list'))
    #"""
    #from py_util.cmdline import run, wrap
    #run(argv, {
        #None: wrap([str, 'safe_eval', 'safe_eval'], {}
            #)( create_dpcofg_table_files),
        #'to_tables': wrap(['eval'], {}
            #)( to_tables),
        #'groupname': wrap([str, 'eval', 'load'], {}
            #)( create_dpcofg_table_files_using_clip_dict),
        #'using_clip_dict': wrap([str, 'eval', 'load'], {}
            #)( create_dpcofg_table_files_using_clip_dict),
        #})


#if __name__ == '__main__':
    #import sys
    #main(sys.argv)


