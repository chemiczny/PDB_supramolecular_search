"""
Created on Sat Apr 21 13:49:15 2018

@author: michal
"""

import numpy as np
from math import sin, cos


def normalize(v):
    """
    Funkcja pomocnicza. Normalizuje wektor v.

    Wejscie:
    v - numpy numeric array, wektor do normalizacji

    Wyjcie:
    v - znormalizowany wektor v
    """
    norm = np.linalg.norm(v)
    if norm == 0:
        return v
    return v / norm


def rotate_vector(vector, axis, angle):
    """
    Obroc wspolrzedne wokol zadanej osi o zadany kat. Czyli wygeneruj macierz obrotu i
    przemnoz wspolrzedne przez nia.
    Wejscie:
    coords - lista wspolrzednych atomow
    norm_vec - os obroty
    angle - kat w radianach!!!

    Wyjscie:
    newCoords - obrocone wspolrzedne
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
