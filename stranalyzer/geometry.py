import numpy 
from numpy.linalg import norm

def delt(xi, xj):
    """
    Compute oriented vector from the first to the second joint
    """
    return xj - xi

def len(xi, xj):
    """
    Compute distance from the first to the second joint
    """
    return  norm(delt(xj, xi))
