import networkx
import sqlite3

def find_subproblems(db):
    g = build_graph(db)
    return nx.connected_components(g)

def build_graph(db):
    g = nx.Graph()
    [g.add_node(x) for (x,) in db.execute("""select id from transcripts""")]
    mids = [x for (x,) in db.execute("""select id from multiplicities""")]
    for mid in mids:
        pairs = db.execute("""select a.transcript,b.transcript 
                              from multiplicity_entries as a
                              join multiplicity_entries as b
                              on a.multiplicity=? and b.multiplicity=?
                              and a.transcript < b.transcript""",
                           (mid,mid))
        [g.add_edge(x,y) for (x,y) in pairs]
    return g
