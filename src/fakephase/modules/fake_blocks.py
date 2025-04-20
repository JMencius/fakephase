import logging
import sys


def dfs(edges: dict, actions: dict, target_site: int) -> int:
    """Depth-First Search (DFS)"""
    for (u, v), w in edges.items():
        if u == target_site:
            if (v not in actions):
                return w * dfs(edges, actions, v)
            else:
                return w * actions[v]

        if v == target_site:
            if (u not in actions):
                return w * dfs(edges, actions, u)
            else:
                return w * actions[u]


def construct_actions(edge: dict, actions: dict, block: list) -> dict:
    addition = None
    while addition != 0:
        for (i, j), a in edge.items():
            addition = 0
            if (i in actions) and (j not in actions):
                actions[j] = actions[i] * a
                addition += 1
            if (j in actions) and (i not in actions):
                actions[i] = actions[j] * a
                addition += 1

    block_variant = set(block)
    cleaned_actions = dict()
    for i in actions:
        if i in block_variant:
            cleaned_actions[i] = actions[i]
    

    return cleaned_actions    


    


def fake_blocks(ref_len: dict, maxratio: float, working_chr: str, variants: dict, triangled_blocks: list, filtered_edges: dict, start_end: tuple) -> dict:
    # find largest two blocks in two near-telomere regions
    s, e = start_end
    regions = [(e - ref_len[working_chr] * maxratio - 1, e)]
    telomere_blocks = list()
    for r in regions:
        start, end = r
        min_length = 2
        for b in triangled_blocks:
            if start <= min(b) < max(b) <= end:
                if len(b) >= min_length:
                    telomere_blocks.append(b)
        
    
    telomere_blocks.sort(key = lambda K : min(K))
    
    # find two closest site to the middle of the chromosome
    min_site = min(variants.keys())

    total_actions = dict()
    for t in [telomere_blocks]:
        count = 0
        for b in t:
            actions = dict()
            # parse action based on edges
            smallest_site = min(b)
        
            # set PS tag to the position of the smallest variant
            if count == 0:
                ps_tag = min_site
                actions[min_site] = 1
            else:
                ps_tag = smallest_site
            actions[smallest_site] = 1
            flag = 1
            ## for actions: 1 means keep the original GT tag, -1 means flip the original GT tag   
        
            actions = construct_actions(filtered_edges, actions, b)
        
            # add ps_tag
            for i in actions:                
                actions[i] = (actions[i], ps_tag)

                 
            total_actions.update(actions)
            count += 1
    
    return total_actions

