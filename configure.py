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
    
    if not isfile(configurationFileName):
        return { "N" : 1, "cif" : "cif/*.cif", "scratch" : "scr" }
        
    configFile = open(configurationFileName)
    config = json.load(configFile)
    configFile.close()
    
    if "scratch" in config:
        if not isdir(config["scratch"]):
            makedirs(config["scratch"])
    else:
        config["scratch"] = "scr"
        
    if "externalLibsPath" in config:
        if not config["externalLibsPath"] in sys.path:
            sys.path.insert(0, config["externalLibsPath"] )
        
    return config
    
    