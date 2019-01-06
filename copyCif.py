from glob import glob
from shutil import copyfile
from os.path import join, basename

cif = glob("/net/archive/groups/plggsuprm/cif/*cif")

for c in cif:
	destiny = join("cif",basename(c))
	copyfile(c, destiny)
