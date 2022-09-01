import networkx as nx
def find_path(graph, start, end, path=[]):
    """ Find path between two nodes in a Graph.

    Parameters
    ----------
    graph: networkx Graph or DiGraph
        Graph representing AA structure of PEI 
    start: networkx node
        Starting node of a path
    end: networkx node
        Last node of a path
    path: list (Optional)
        List containing nodes in the path. (default [])

    Returns
    -------
???
    """
    path = path + [start] # for path=[], path is now [start]
    if start==end: #Zero length path.
        return path
    if not graph.has_key(start): #No path.
        return None
    for node in graph[start]:
        if node not in path:
            newpath=find_path(graph,node,end,path)
            if newpath: 
                return newpath
    return None


G=nx.Graph()
G.add_edges_from([(0,1),(1,2),(2,3),(0,4),(4,3)])
path=find_path(G,0,3)
print(path)
