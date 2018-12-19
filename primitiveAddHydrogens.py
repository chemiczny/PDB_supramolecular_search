import sys
from pymol import cmd

pdbFilename = sys.argv[1]
cmd.load(pdbFilename)
cmd.h_add("(all)")
cmd.save(pdbFilename, "(all)")