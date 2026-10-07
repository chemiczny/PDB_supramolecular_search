# PDB_supramolecular_search

Tools for finding and analysing interactions involving aromatic rings (anion-pi, cation-pi, methyl-pi, pi-pi) in PDB structures stored as mmCIF files. 

## Repository layout
```
supramolecular_search/   core library (CIF analysis, ring detection, anion recognition, logging)
  anion_templates/       anion templates used by the anion recogniser
scripts/                 command-line entry points (simple_run.py, brutal_*.py, ...)
hpc/                     example SLURM job files for the scripts
pymol_plugin/            PyMOL plugin for browsing and filtering results
benchmarks/              performance comparisons
tests/                   pytest suite
```

## Installation
Python 3.10+ is required. Install the package (editable mode is convenient for development):

```
pip install -e .
```

Optional extras: `pip install -e ".[analysis]"` (matplotlib, pandas, requests for the post-processing scripts)
and `pip install -e ".[test]"` (pytest).

The PyMOL plugin has its own requirements in pymol_plugin/requirements.txt. Instructions for installing PyMOL plugins can be found here: https://pymolwiki.org/index.php/Plugins

## Getting started
All scripts read `config.json` and write their results relative to the current working directory,
so run them from the directory that holds your configuration (e.g. the repository root).
The easiest way to create the configuration is to copy the example file and modify it:

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
These 3 parameters indicate: how many processes are used for parsing, where to read the CIF files from and where to store temporary files.
The simplest way to run the analysis is to use simple_run.py:
```
python3 scripts/simple_run.py
```
Then the result files should be available in the logs directory:
```
ls logs
additionalInfo.log  anionPi.log   cif2process.log  linearAnionPi.log  methylPi.log  planarAnionPi.log    timeStart.log
anionCation.log     cationPi.log  hBonds.log       metalLigand.log    piPi.log      progressSummary.log  timeStop.log
```
If you have a large number of CIF files (for example the entire PDB), you will probably want to run the analysis on
a cluster. For that, use scripts/brutal_prepare.py, scripts/brutal_run.py and scripts/brutal_finish.py, in this order
(brutal_prepare.py must be run every time before brutal_run.py). Example SLURM files are in the hpc directory; submit
them from the directory containing config.json, e.g. `sbatch hpc/brutal_prepare.slurm`.

## Running tests
```
pip install -e ".[test]"
pytest
```

## Pymol plugin
After installation, you should be able to see "Supramolecular analyser" among other plugins. Well, it is written with Tkinter, so it looks like software from 
90s, but it works even better.

Anyway, it is possible to load any result log independently for each tab, or to read entire directory with logs (with original file names like anionPi.log etc).
After that filtering, sorting, merging, excluding results can be done with one click. Interface is a little bit unintuitive. In short: always pay attention to selected checkbox.
