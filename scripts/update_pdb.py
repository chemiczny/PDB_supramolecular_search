from supramolecular_search.config import configure

config = configure()

from Bio.PDB import PDBList
from os.path import dirname, join, basename
from glob import glob
from os import remove

pdb_dir = dirname(config["cif"])
pl = PDBList(pdb=pdb_dir)
pl.flat_tree = True

existing_cifs = glob(join(pdb_dir, "*.cif"))
existing_pdbs = set([])
for cif in existing_cifs:
    existing_pdbs.add(basename(cif).replace(".cif", ""))

all_pdbs = set(pdb_code.lower() for pdb_code in pl.get_all_entries())

pdb2delete = existing_pdbs - all_pdbs
pdb2download = all_pdbs - existing_pdbs

print("Found ", len(pdb2delete), " files to delete")
print("and ", len(pdb2download), " file to download")

for pdb_code in pdb2download:
    pl.retrieve_pdb_file(pdb_code, file_format="mmCif")

for pdb_code in pdb2delete:
    cif_path = join(pdb_dir, pdb_code + ".cif")
    remove(cif_path)
