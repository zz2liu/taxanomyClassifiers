#!/usr/bin/env python
"""get the rdp classification detailed data from RDP classifier web tool.
"""
from time import sleep
from pdb import set_trace
from py_util.cbook.cookie_url import CookieOpener
from shutil import copyfileobj
from numpy import disp
from time import time


rdp_base = 'http://rdp.cme.msu.edu/classifier/'
servlet = 'classifierServlet' #the form action url
check_page = 'cl_status.jsp' #the refresh url

def get_download_url(lines):
    """gets download url (relative) from lines"""
    target = 'detail.jsp?root='
    #lines = list(lines) #debug
    for line in lines:
        if target in line:
            fields = line.split("'", 2)
            root_url = fields[1]
            return root_url.replace('detail', 'download')
    else:
        set_trace()
        raise ValueError('%s not found.' % target)


def is_processing(result):
    """result is consumed to the title line."""
    for line in result:
        if '<title>' in line:
            if 'Classifier Status</title>' in line:
                return True
            else: #break here
                return False

def download_assignments(opener, fasta_fname, interval=3):
    """download the rdp assignments to each seq in fasta.
    
    - interval=3: the interval between each refresh.
    """
    params = {"file" : open(fasta_fname, "rb") }
    #submit and refresh until processed
    result = opener.open(rdp_base+servlet, params)
    while is_processing(result):
        sleep(interval)
        result = opener.open(rdp_base + check_page)

    #download the detailed text result
    result = opener.open(rdp_base + get_download_url(result))
    return result

def classify_fastas(fastas, out_fnames=None, interval=3):
    """classify each fasta and save result to each out_fname.

    - fastas: a list of fasta file names.
    - out_fnames=None: a list of out_filenames for output.
    - interval=3: the refresh interval in second.
    """
    opener = CookieOpener()
    if not out_fnames:
        out_fnames = ['%s.rdp' % name for name in fastas]
    for in_fname, out_fname in zip(fastas, out_fnames):
        result = download_assignments(opener, in_fname, interval)
        copyfileobj(result, open(out_fname, 'wb'))
        disp('%s --> %s' % (in_fname, out_fname))

def timetest(fastas, out_fnames=None, interval=1):
    """classify each fasta and save result to each out_fname.

    - fastas: a list of fasta file names.
    - out_fnames: a list of out_filenames for output.
    - interval: the refresh interval in second.
    """
    opener = CookieOpener()
    if not out_fnames:
        out_fnames = ['%s.rdp' % name for name in fastas]
    result = []
    for in_fname, out_fname in zip(fastas, out_fnames):
        start = time()
        curr = download_assignments(opener, in_fname, interval)
        copyfileobj(curr, open(out_fname, 'wb'))
        result.append(time()-start)
        disp('%s --> %s' % (in_fname, out_fname))
    return result

def main(argv):
    """usage: prog.py <a list of fasta names> [--interval=int].

    eg: p250 = Popen('python get_rdp.py \"%s\"' %
    glob('clip_250/*.fa'), stdout=open('rdp250.log', 'w'),
    shell=True)
    """
    from py_util.cmdline import run, wrap
    run(argv, {
        None: wrap(['safe_eval', 'safe_eval', int],
            {'interval': int})(
            classify_fastas),
        })
    return

if __name__ == '__main__':
    import sys; main(sys.argv)
    #filename = '/home/zongzhi/Working/TaxonomyFinishing/test.fa'

