class UnionFind:
    """Union-Find data structure for merging nodes that satisfy triangular consistency."""
    
    def __init__(self):
        self.parent = {}

    def find(self, node):
        """Find the representative of the set with path compression."""
        if self.parent[node] != node:
            self.parent[node] = self.find(self.parent[node])
        return self.parent[node]

    def union(self, node1, node2):
        """Merge two nodes into the same set."""
        root1 = self.find(node1)
        root2 = self.find(node2)
        if root1 != root2:
            # Merge into the same set
            self.parent[root2] = root1

    def add(self, node):
        """Initialize a node in the Union-Find structure."""
        if node not in self.parent:
            self.parent[node] = node

    def get_connected_components(self):
        """Retrieve all connected components."""
        groups = {}
        for node in self.parent:
            root = self.find(node)
            if root not in groups:
                groups[root] = set()
            groups[root].add(node)
        return list(groups.values())
