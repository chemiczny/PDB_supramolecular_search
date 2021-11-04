# PDB_supramolecular_search

This repository contains a bunch of scripts for analyzing interactions involving aromatic ring (anion-pi, cation-pi, methyl-pi, pi-pi) in cif files. 

## Requirements
All required external modules are listed in requirements.txt (for scripts related to analyzing cif files) and pymol_plugin/requirements.txt (for pymol plugin).

## Installation
Just download files, and use it. Instructions for installing pymol plugin can be found here: https://pymolwiki.org/index.php/Plugins

## Getting started
Firstly, you have to create configuration file. The easiest way is just to copy example configuration file and then modify it.

```
cp configExample.json config.json
```

This example config.json looks like this:
```json
{
	"N" : 6,
	"cif" : "cif/*.cif",
	"scratch" : "scratch"
}
```
These 3 parameters indicates: how many processors are selected for parsing, where read the files from and where store the temporary files.
The simplest way to run analysis is to use simpleRun.py.
```
python3 simpleRun.py
```
Then result files should be available in logs directory.
```
ls logs
additionalInfo.log  anionPi.log   cif2process.log  linearAnionPi.log  methylPi.log  planarAnionPi.log    timeStart.log
anionCation.log     cationPi.log  hBonds.log       metalLigand.log    piPi.log      progressSummary.log  timeStop.log
```
If you have large amount of cif files (for example entire PDB database), you may probably want to run analysis on some kind
of supercomputer. For such a situation, there are other scripts: brutalPrepare.py, brutalRun.py and brutalFinish.py. Example slurm files 
executing them are also available: brutalPrepare.slurm, brutalRun.slurm and brutalFinish.slurm. These scripts must be run in metioned order (brutalPrepare.py must be run
everytime before brutalRun.py).

##Pymol plugin
After installation, you should be able to see "Supramolecular analyser" among other plugins. Well, it is written with Tkinter, so it looks like software from 
90s, but it works even better.

Anyway, it is possible to load any result log independently for each tab, or to read entire directory with logs (with original file names like anionPi.log etc).
After that filtering, sorting, merging, excluding results can be done with one click. Interface is a little bit unintuitive. In short: always pay attention to selected checkbox.
