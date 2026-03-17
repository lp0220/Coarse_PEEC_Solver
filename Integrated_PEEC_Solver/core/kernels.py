import numpy as np
import math
from config import E0, U0


def cal_p_self(b):
    """Calculate self coefficient of potential for a node."""
    return 2 * math.log(1 + math.sqrt(2)) / (math.pi * E0 * b)


def cal_p_oth(p1, p2, b):
    """Calculate mutual coefficient of potential between two nodes."""
    dist = np.linalg.norm(p1 - p2)
    return 1 / (dist * 4 * math.pi * E0)


def cal_l_oth_approx(points, c1, c2):
    """
    Calculate mutual inductance between two branches (approximation).
    c1, c2: branch node pairs, 1-based indices.
    """
    pa1, pa2 = points[c1[0] - 1], points[c1[1] - 1]
    pb1, pb2 = points[c2[0] - 1], points[c2[1] - 1]

    mid1 = (pa1 + pa2) / 2
    mid2 = (pb1 + pb2) / 2
    dist = np.linalg.norm(mid1 - mid2)

    v1 = (pa2 - pa1)
    v2 = (pb2 - pb1)
    return U0 * np.dot(v1, v2) / (dist * 4 * math.pi)
