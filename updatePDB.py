from configure import configure
config = configure()

from Bio.PDB import PDBList
from os.path import dirname, join, basename
from glob import glob
from os import remove

pdb_dir = dirname(config["cif"])
pl = PDBList(pdb=pdb_dir)
pl.flat_tree = True

existingCifs = glob( join( pdb_dir, "*.cif" ) )
existingPdbs = set( [  ] )
for cif in existingCifs:
	existingPdbs.add( basename(cif).replace(".cif", "") )

allPdbs = set( pdb_code.lower() for pdb_code in  pl.get_all_entries())

pdb2delete = existingPdbs - allPdbs
pdb2download = allPdbs - existingPdbs

print("Found ", len(pdb2delete), " files to delete")
print("and ", len(pdb2download), " file to download")

for pdb_code in pdb2download: 
	pl.retrieve_pdb_file(pdb_code, file_format="mmCif") 

for pdb_code in pdb2delete:
	cifPath = join(pdb_dir, pdb_code+".cif")
	remove(cifPath)