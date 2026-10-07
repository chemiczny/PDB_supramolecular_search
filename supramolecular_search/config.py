#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 22:46:10 2020

@author: michal
"""
import json
from os.path import isfile, isdir
from os import makedirs
import sys

def configure():
    configurationFileName = "config.json"
    
    if isfile(configurationFileName):
        configFile = open(configurationFileName)
        config = json.load(configFile)
        configFile.close()
    else:
        config = { "N" : 1, "cif" : "cif/*.cif" }

    if not "scratch" in config:
        config["scratch"] = "scr"

    if not isdir(config["scratch"]):
        makedirs(config["scratch"])
        
    if "externalLibsPath" in config:
        if not config["externalLibsPath"] in sys.path:
            sys.path.insert(0, config["externalLibsPath"] )
        
    return config
    
    