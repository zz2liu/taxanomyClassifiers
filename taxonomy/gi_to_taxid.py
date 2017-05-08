"""gi_to_taxid.py

06/28/07 functions to a class.

todo: error handler for .connect
"""
from pdb import set_trace
import sqlite3

table_name = 'taxids' #the main table in db.
def create_db(db_path, gi_taxid_dmp):
    """ 
    - db_path:
    - gi_taxid_dmp: path of .dmp file from ncbi
    """
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    cur.execute("""create table %s
            (gi integer primary key, taxid integer)
            """ % table_name)
    INSERT = "insert into %s values (?, ?)" % table_name
    for i, line in enumerate(file(gi_taxid_dmp)):
        if not (i % 100000): print '.',
        cur.execute(INSERT, map(int, line.split()))
    con.commit()
    con.close()

class TaxidFinder(object):
    """

    Note: you must create a (gi, taxid) database before your initiate.
    """
    createDb = staticmethod(create_db)

    def __init__(self, db_path):
        """
        - db_path: the exist gi_taxid.db with fields (gi integer, taxid integer)
        """
        self._con = sqlite3.connect(db_path)
        self._cur = self._con.cursor()

    def __del__(self):
        self._con.close()


    def find(self, gis, return_more=False):
        """
        - gis: a seq of ncbi gi.
        - return_more=False: if true, also return queries not found.
        """
        cur = self._cur
        ## gis into temp table
        cur.execute("""create temp table gis (gi_ integer) """)
        for gi in gis:
            cur.execute("""insert into gis values (?)""", (gi,))

        cur.execute("""select gi, taxid from taxids
                inner join gis on gi==gi_
                """)
        result = [row for row in cur]

        if return_more:
            cur.execute("""select gi_ from gis
                where not exists (
                    select * from taxids
                    where taxids.gi==gis.gi_)
                """)
            bad_gis = [row[0] for row in cur]

        cur.execute("""drop table gis""")
        if return_more:
            return result, bad_gis
        else:
            return result

    __call__ = find

