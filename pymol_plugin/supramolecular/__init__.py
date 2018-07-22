#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 22 19:35:22 2018

@author: michal
"""
import sys
import os
path = os.path.dirname(__file__)
if not path in sys.path:
    sys.path.append(path)
    
from fetchDialog import fetchdialog

if sys.version_info[0] < 3:
    from pymol import plugins
else:
    pass
    
def __init_plugin__(self=None):
    plugins.addmenuitem('Supramolecular analyser', fetchdialog)
   
if __name__ == "__main__":
    fetchdialog(True)