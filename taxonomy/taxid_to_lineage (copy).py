"""taxid_to_lineage.py
"""
import sqlite3

def create_db(nodes_dmp, names_dmp, db_path):
    """create the db with two tables (nodes and names).

    - nodes_dmp, names_dmp: lines.
    """
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    ## init table nodes
    cur.execute("""create table nodes
            (taxid integer primary key, parent integer, rank text)
            """)
    INSERT = "insert into nodes values (?, ?, ?)"
    for line in nodes_dmp:
        taxid, parent, rank = line.split('\t|\t')[:3]
        cur.execute(INSERT, (int(taxid), int(parent), rank))

    ## init table names as [(taxid, scientific_name)]
    cur.execute("""create table names
            (taxid integer primary key, name text)
            """)
    INSERT = "insert into names values (?, ?)"
    for line in names_dmp:
        taxid, name, unique_name, name_class = line.split('\t|\t')
        if name_class.startswith('scientific name'):
            cur.execute(INSERT, (int(taxid), name))
    con.commit()
    con.close()



def get_id_rank_names(cur, taxids):
    cur.execute("""create temp table temp (taxid integer primary key)""")
    for id in taxids:
        cur.execute("""insert into temp values (?)""", (id,))

    cur.execute("""select * from temp""")
    cur.execute("""select temp.taxid, nodes.rank, names.name
            from temp
            inner join nodes on nodes.taxid==temp.taxid
            inner join names on names.taxid==temp.taxid
            """)
    result = [row for row in cur]
    cur.execute("""drop table temp""")
    return result

def get_ancestor_ids(cur, node_id):
    """return [taxid] of self and ancestors."""
    result = [node_id]
    while node_id != 1:
        parent = _get_parent_id(cur, node_id)
        result.append(parent)
        node_id = parent
    return result

def _get_parent_id(cur, node_id):
    """return parent_id."""
    cur.execute("""select parent from nodes
            where taxid == ?
            """, (node_id,))
    return cur.next()[0]


def taxid_to_lineage(cur, taxid):
    """return a [(rank, name)] in leaf to root order for taxid.
    """
    ancestor_ids = get_ancestor_ids(cur, taxid)
    return get_id_rank_names(cur, ancestor_ids)

def taxids_to_lineages(taxids, db_path):
    """return a [(rank, name)] in leaf to root order for each taxid.
    """
    con = sqlite3.connect(db_path)
    cur = con.cursor()
    for id in taxids:
        yield id, taxid_to_lineage(cur, id)
    con.close()
