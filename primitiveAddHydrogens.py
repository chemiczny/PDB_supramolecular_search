import sys
from pymol import cmd

pdbFilename = sys.argv[1]
cmd.load(pdbFilename)
cmd.h_add("(element O or element N)")
cmd.save(pdbFilename, "(all)")
