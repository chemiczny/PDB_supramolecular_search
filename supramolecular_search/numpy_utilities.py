#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 13:49:15 2018

@author: michal
"""

import numpy as np
from math import sin, cos, radians

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

def rotateVector( vector, axis, angle ):
    """
    Obroc wspolrzedne wokol zadanej osi o zadany kat. Czyli wygeneruj macierz obrotu i
    przemnoz wspolrzedne przez nia. 
    Wejscie:
    coords - lista wspolrzednych atomow
    normVec - os obroty
    angle - kat w radianach!!!
    
    Wyjscie:
    newCoords - obrocone wspolrzedne
    """
    normVec = normalize(axis)
    a,b,c = normVec
    cosDiff = 1-cos(angle)
    cosA = cos(angle)
#    print(angle*180./3.1416)
    rotateMatrix = np.array( [ [ cosA+ a*a*cosDiff, a*b*cosDiff - c * sin(angle), a*c*cosDiff + b *sin(angle)   ] ,
                            [a*b*cosDiff + c * sin(angle),cosA+ b*b*cosDiff,b*c*cosDiff - a * sin(angle)],
                            [a*c*cosDiff - b * sin(angle),c*b*cosDiff + a * sin(angle),cosA+ c*c*cosDiff ]  ] )

    return rotateMatrix.dot(vector)