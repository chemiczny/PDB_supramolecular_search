#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 13:49:15 2018

@author: michal
"""

import numpy as np

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