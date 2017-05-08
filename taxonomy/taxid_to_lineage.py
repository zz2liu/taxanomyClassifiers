"""taxid_to_lineage.py  -- using sqlite

todo: prepare get_parent
"""
import os
from warnings import warn
import sqlite3

def create_db(nodes_dmp, names_dmp, db_path):
    """create the db with two tables (nodes and names).

    - nodes_dmp, names_dmp: lines or filenames.
    """
    if isinstance(nodes_dmp, str):
        nodes_dmp = file(nodes_dmp)
    if isinstance(names_dmp, str):
        names_dmp = file(names_dmp)
    if os.path.isfile(db_path):
        os.remove(db_path)

    con = sqlite3.connect(db_path)
    cur = con.cursor()
    ## init table nodes
    cur.execute("""create table nodes
            (taxid integer primary key, parent integer, rank text)
            """)
    INSERT = "insert into nodes values (?, ?, ?)"
    for line in nodes_dmp:
        try:
            taxid, parent, rank = line.split('\t|\t')[:3]
        except ValueError, e:
            warn("%s\n%s" % (line, e))
            continue
        cur.execute(INSERT, (int(taxid), int(parent), rank))

    ## init table names as [(taxid, scientific_name)]
    cur.execute("""create table names
            (taxid integer primary key, name text)
            """)
    INSERT = "insert into names values (?, ?)"
    for line in names_dmp:
        try:
            taxid, name, unique_name, name_class = line.split('\t|\t')
        except ValueError, e:
            warn("%s\n%s" % (line, e))
            continue
        if name_class.startswith('scientific name'):
            cur.execute(INSERT, (int(taxid), name))
    con.commit()
    con.close()


class LineageFinder(object):
    """
    Note: you should built a valid database, before initiate.
    """
    createDb = staticmethod(create_db)

    def __init__(self, db_path):
        """
        - db_path: the path to the sqlite3 db containing tables of names and
          nodes.
        """
        self._con = sqlite3.connect(db_path)
        self._cur = self._con.cursor()

    def __del__(self):
        self._con.close()

    def findLineages(self, taxids):
        """return a [(rank, name)] in leaf to root order for each taxid.

        - taxids: a seq of taxid as integers.
        """
        cur = self._cur
        for id in taxids:
            yield id, self._taxid_to_lineage(cur, id)
    __call__ = findLineages

    @staticmethod
    def _taxid_to_lineage(cur, taxid):
        """return a [(rank, name)] in leaf to root order for taxid.
        """
        ancestor_ids = LineageFinder._get_ancestor_ids(cur, taxid)
        return LineageFinder._get_id_rank_names(cur, ancestor_ids)

    @staticmethod
    def _get_id_rank_names(cur, taxids):
        cur.execute("""create temp table temp (taxid integer primary key)""")
        for id in taxids:
            cur.execute("""insert into temp values (?)""", (id,))

        cur.execute("""select * from temp""")
        cur.execute("""select temp.taxid, nodes.rank, names.name
                from temp
                inner join nodes on nodes.taxid==temp.taxid
                inner join names on names.taxid==temp.taxid
                """)
        result = [row[1:] for row in cur]
        cur.execute("""drop table temp""")
        return result

    @staticmethod
    def _get_ancestor_ids(cur, node_id):
        """return [taxid] of node and ancestors."""
        result = [node_id]
        while node_id != 1:
            parent = LineageFinder._get_parent_id(cur, node_id)
            result.append(parent)
            node_id = parent
        return result

    @staticmethod
    def _get_parent_id(cur, node_id):
        """return parent_id."""
        cur.execute("""select parent from nodes
                where taxid == ?
                """, (node_id,))
        return cur.next()[0]


