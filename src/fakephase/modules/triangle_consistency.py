import logging


def triangle_consistency(sign_edges: dict, mincoverage: int, conf: float) -> dict:
    # prefilter connections based on coverage and confidence
    filtered_sign_edges = dict()
    for k, j in sign_edges.items():
         if (j[0] + j[1]) > mincoverage:
            if max(j[0], j[1]) / (j[0] + j[1]) > conf:
                if j[0] > j[1]:
                    filtered_sign_edges[k] = 1
                if j[0] < j[1]:
                    filtered_sign_edges[k] = -1

    # triangle consistency check
    fit_sign_edges = dict()
    for k, j in filtered_sign_edges.items():
        
