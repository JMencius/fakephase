import logging
from fakephase.classes.unionfind import UnionFind


def build_phase_blocks(sign_edges: dict, mincoverage: int, conf: float) -> dict:
    # prefilter connections based on coverage and confidence
    filtered_sign_edges = dict()
    for k, j in sign_edges.items():
         if (j[0] + j[1]) > mincoverage:
            if max(j[0], j[1]) / (j[0] + j[1]) > conf:
                if j[0] > j[1]:
                    filtered_sign_edges[k] = 1
                if j[0] < j[1]:
                    filtered_sign_edges[k] = -1

    # generate adjacent dict
    adj_dict = dict()
    for u, v in filtered_sign_edges.keys():
        adj_dict.setdefault(u, set()).add(v)
        adj_dict.setdefault(v, set()).add(u)

    uf = UnionFind()
    consistent_triangles = set()

    # loop every possible triangle
    for a in adj_dict:
        neighbors = list(adj_dict[a])
        n = len(neighbors)
        for i in range(n - 1):
            for j in range(i + 1, n):
                b, c = neighbors[i], neighbors[j]
                if tuple(sorted([b, c])) in edge_set:
                    triangle = tuple(sorted([a, b, c]))
                    if is_triangle_consistent(triangle, filtered_sign_edges):
                        consistent_triangles.add(triangle)
                        for node in triangle:
                            uf.add(node)
                        uf.union(a, b)
                        uf.union(a, c)
                        uf.union(b, c)                   
    
    blocks = uf.get_connected_components()
    return blocks


def is_triangle_consistent(triangle: tuple, edge_dict: dict) -> bool:
    pair1 = tuple(sorted([triangle[0], triangle[1]]))
    pair2 = tuple(sorted([triangle[0], triangle[2]]))
    pair3 = tuple(sorted([triangle[1], triangle[2]]))
    
    if edge_dict[pair1] * edge_dict[pair2] == edge_dict[pair3]:
        return True
    else:
        return False                  
