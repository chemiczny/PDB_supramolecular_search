#!/usr/bin/env python3
"""
Created on Mon Oct  1 14:26:23 2018

@author: michal
"""

from supramolecular_search.supramolecular_logging import (
    write_anion_pi_header,
    write_anion_cation_header,
    write_pi_pi_header,
    write_cation_pi_header,
    write_hbonds_header,
    write_metal_ligand_header,
)
from supramolecular_search.supramolecular_logging import (
    write_anion_pi_linear_header,
    write_anion_pi_planar_header,
    write_methyl_pi_header,
)
from os.path import isdir, join
from os import makedirs, remove
import glob
import time
from supramolecular_search.config import configure

config = configure()

cif_pattern = config["cif"]

if not isdir("logs"):
    makedirs("logs")

write_anion_pi_header()
write_anion_cation_header()
write_pi_pi_header()
write_cation_pi_header()
write_hbonds_header()
write_metal_ligand_header()
write_anion_pi_linear_header()
write_anion_pi_planar_header()
write_methyl_pi_header()

open("logs/additionalInfo.log", "w").close()

time_start = time.time()
time_file = open("logs/timeStart.log", "w")
time_file.write(str(time_start))
time_file.close()


cif_files = glob.glob(cif_pattern)
cif_no_file = open("logs/cif2process.log", "w")
cif_no_file.write(str(len(cif_files)))
cif_no_file.close()

scratch = config["scratch"]
files2remove = glob.glob(join(scratch, "*"))

for f in files2remove:
    remove(f)

files_for_step = 70
actual_id = 0
files_for_current_step = 0

input_file = open("cif2process.dat", "w")

actual_file = open(join(scratch, str(actual_id) + ".dat"), "w")
input_file.write(join(scratch, str(actual_id) + ".dat") + "\n")
for cif in cif_files:
    if files_for_current_step > files_for_step:
        actual_id += 1
        actual_file.close()
        actual_file = open(join(scratch, str(actual_id) + ".dat"), "w")
        input_file.write(join(scratch, str(actual_id) + ".dat") + "\n")
        files_for_current_step = 0

    actual_file.write(cif + "\n")
    files_for_current_step += 1

actual_file.close()
input_file.close()
