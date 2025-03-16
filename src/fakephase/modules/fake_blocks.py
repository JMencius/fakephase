import loggings
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


def find_closest_variants(variant_sites: list, mid_of_chrom: int) -> tuple:
    dist = list()
    for i in variant_sites:
        dist.append((i, abs(i - mid_of_chrom)))
    dist.sort(key = lambda K : K[1])
    
    return (dist[0][0], dist[1][0])
    


def fake_blocks(ref_len: dict, working_chr: str, variants: dict, triangled_blocks: list, edges: dict, max_len: int) -> dict:
    # find largest two blocks in two near-telomere regions
    regions = [(0, maxlen + 1), (ref_len[working_chr] - maxlen - 1, ref_len[working_chr])]
    longest_blocks = list()
    for r in regions:
        start, end = r
        length = 0
        for b in triangled_blocks:
            if start <= min(b) < max(b) <= end:
                if len(b) > length:
                    length = len(b)
                    max_block = b
        if length == 0:
            logging.error("Can not find input block for fake_blocks")
            sys.exit(1)
        else:
            longest_blocks.append(max_block)
    
    longest_blocks.sort(key = lambda K : min(K))
    
    # find two closest site to the middle of the chromosome
    mid_sites = find_closest_variants(variants.keys(), ref_len[working_chr] // 2)

    actions = dict()
    c = 0
    for b in longest_blocks:
        # parse action based on edges
        smallest_site = min(b)
        
        # set PS tag to the position of the smallest variant
        ps_tag = smallest_site
        fake_blocks.append((actions, ps_tag))

        ## for actions: 1 means keep the original GT tag, -1 means flip the original GT tag   
        actions[smallest_site] = 1
        for site in b:
            if site not in actions:
                actions[site] = (dfs(edges, actions, site), ps_tag)                      

        actions[mid_sites[c]] = (1, ps_tag)


        # index add 1 for mid_sites
        c += 1

    return actions
