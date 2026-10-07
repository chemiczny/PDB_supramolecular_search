"""
Created on Thu Apr  2 22:46:10 2020

@author: michal
"""

import json
from os.path import isfile, isdir
from os import makedirs
import sys


def configure():
    configuration_file_name = "config.json"

    if isfile(configuration_file_name):
        config_file = open(configuration_file_name)
        config = json.load(config_file)
        config_file.close()
    else:
        config = {"N": 1, "cif": "cif/*.cif"}

    if "scratch" not in config:
        config["scratch"] = "scr"

    if not isdir(config["scratch"]):
        makedirs(config["scratch"])

    if "externalLibsPath" in config:
        if config["externalLibsPath"] not in sys.path:
            sys.path.insert(0, config["externalLibsPath"])

    return config
