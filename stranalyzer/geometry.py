from numpy import diff

def delt(xi, xj):
    """
    Compute oriented vector from the first to the second joint
    """
    return xj - xi
