"""Various utility functions"""

def get_oldest_mcp_parent(mcp, n_iters=0):
    """Recursively looks for the oldest parent of the input MCParticle"""
    parents = mcp.getParents()
    if (len(parents) < 1):
        return mcp, n_iters
    for parent in parents:
        # Skipping if the particle is its own parent
        if parent is mcp:
            continue
        return get_oldest_mcp_parent(parent, n_iters+1)
