"""
Created on Sat Apr 21 13:49:15 2018

@author: michal
"""

import numpy as np
from math import sin, cos


def normalize(v):
    """Return v scaled to unit length.

    Args:
        v: numeric numpy array.

    Returns:
        The normalized vector.
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def rotate_vector(vector, axis, angle):
    """Rotate a vector around an axis by the given angle.

    Args:
        vector: 3-element vector to rotate.
        axis: rotation axis (does not need to be normalized).
        angle: rotation angle in radians.

    Returns:
        The rotated vector.
    """
    norm_vec = normalize(axis)
    a, b, c = norm_vec
    cos_diff = 1 - cos(angle)
    cos_a = cos(angle)
    rotate_matrix = np.array(
        [
            [
                cos_a + a * a * cos_diff,
                a * b * cos_diff - c * sin(angle),
                a * c * cos_diff + b * sin(angle),
            ],
            [
                a * b * cos_diff + c * sin(angle),
                cos_a + b * b * cos_diff,
                b * c * cos_diff - a * sin(angle),
            ],
            [
                a * c * cos_diff - b * sin(angle),
                c * b * cos_diff + a * sin(angle),
                cos_a + c * c * cos_diff,
            ],
        ]
    )

    return rotate_matrix.dot(vector)
