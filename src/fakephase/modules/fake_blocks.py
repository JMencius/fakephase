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



def find_closest_variants(variant_sites: list, mid_of_chrom: int) -> tuple:
    dist = list()
    for i in variant_sites:
        dist.append((i, abs(i - mid_of_chrom)))
    dist.sort(key = lambda K : K[1])
    
    return (dist[0][0], dist[1][0])
    


def fake_blocks(ref_len: dict, working_chr: str, variants: dict, triangled_blocks: list, filtered_edges: dict, maxlen: int, start_end: tuple) -> dict:
    # find largest two blocks in two near-telomere regions
    s, e = start_end
    regions = [(s, s + maxlen + 1), (e - maxlen - 1, e)]
    telomere_blocks1 = list()
    telomere_blocks2 = list()
    c = 0
    for r in regions:
        start, end = r
        min_length = 2
        for b in triangled_blocks:
            if start <= min(b) < max(b) <= end:
                if len(b) >= min_length:
                    if c == 0:
                        telomere_blocks1.append(b)
                    else:
                        telomere_blocks2.append(b)
        c += 1
        
    
    telomere_blocks1.sort(key = lambda K : min(K))
    telomere_blocks2.sort(key = lambda K : min(K))
    
    # find two closest site to the middle of the chromosome
    mid_sites = find_closest_variants(variants.keys(), ref_len[working_chr] // 2)

    total_actions = dict()
    c = 0
    for t in [telomere_blocks1, telomere_blocks2]:
        flag = 0
        for b in t:
            actions = dict()
            # parse action based on edges
            smallest_site = min(b)
        
            # set PS tag to the position of the smallest variant
            ps_tag = smallest_site
            flag = 1
            ## for actions: 1 means keep the original GT tag, -1 means flip the original GT tag   
            actions[smallest_site] = 1
        
            actions = construct_actions(filtered_edges, actions, b)
        
            # add ps_tag
            for i in actions:                
                actions[i] = (actions[i], ps_tag)                      
            total_actions.update(actions)
        if flag:
            total_actions[mid_sites[c]] = (1, ps_tag)

        # index add 1 for mid_sites
        c += 1
    
    return total_actions
