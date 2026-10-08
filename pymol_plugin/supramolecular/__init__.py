"""
Created on Sun Jul 22 19:35:22 2018

@author: michal
"""

import sys
import os

path = os.path.dirname(__file__)
if path not in sys.path:
    sys.path.append(path)

from fetch_dialog import fetch_dialog

try:
    from pymol import plugins
except ImportError:
    pass


def __init_plugin__(self=None):
    plugins.addmenuitem("Supramolecular analyser", fetch_dialog)


if __name__ == "__main__":
    fetch_dialog(True)
