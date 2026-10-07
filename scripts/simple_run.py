"""
Created on Mon Jan  1 18:00:25 2018

@author: michal
"""

from supramolecular_search.cif_analyser import find_supramolecular
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
from os.path import isdir, basename, join
from os import makedirs, remove
from supramolecular_search.config import configure
import glob
import time
from multiprocessing import Pool
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(processName)s %(levelname)s %(name)s: %(message)s",
)
# import cProfile, pstats, io

###profiler start
# pr = cProfile.Profile()
# pr.enable()

##############

config = configure()
number_of_processes = config["N"]
cif_pattern = config["cif"]
scratch = config["scratch"]

files2remove = glob.glob(join(scratch, "*"))
for f in files2remove:
    remove(f)

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

cif_files = glob.glob(cif_pattern)

data_processed = 0
structure_saved = 0
data_len = len(cif_files)

pool = Pool(number_of_processes)


def prepare_arguments_list(cif_files):
    arguments = []

    for cif in cif_files:
        pdb_code = basename(cif).split(".")[0].upper()
        arguments.append((cif, pdb_code, "default"))

    return arguments


time_start = time.time()
time_file = open("logs/timeStart.log", "w")
time_file.write(str(time_start))
time_file.close()

cif_no_file = open("logs/cif2process.log", "w")
cif_no_file.write(str(len(cif_files)))
cif_no_file.close()


arguments_list = prepare_arguments_list(cif_files)
pool.map(find_supramolecular, arguments_list)
#
# for arg in argumentsList:
#    findSupramolecular(arg)


def merge_logs(log_final, logs):
    final_log = open(log_final, "a+")
    log_files = glob.glob(logs)
    for log_file in log_files:
        if log_file == log_final:
            continue

        new_log = open(log_file, "r")

        line = new_log.readline()
        while line:
            final_log.write(line)
            line = new_log.readline()

        new_log.close()
        remove(log_file)

    final_log.close()


#
def merge_progress_files():
    cif_processed = 0
    log_files = glob.glob(join(scratch, "partialProgress*"))
    for log_file in log_files:
        log = open(log_file, "r")
        cif_processed += int(log.readline())
        log.close()
        remove(log_file)

    cif_no_file = open("logs/cif2process.log", "r")
    cif_no = int(cif_no_file.readline())
    cif_no_file.close()

    progress_summary = open("logs/progressSummary.log", "w")
    progress_summary.write(
        "Przetworzono: " + str(cif_processed) + "/" + str(cif_no) + "\n"
    )
    progress_summary.close()


merge_logs("logs/anionPi.log", join(scratch, "anionPi*.log"))
merge_logs("logs/cationPi.log", join(scratch, "cationPi*.log"))
merge_logs("logs/piPi.log", join(scratch, "piPi*.log"))
merge_logs("logs/anionCation.log", join(scratch, "anionCation*.log"))
merge_logs("logs/hBonds.log", join(scratch, "hBonds*.log"))
merge_logs("logs/metalLigand.log", join(scratch, "metalLigand*.log"))
merge_logs("logs/linearAnionPi.log", join(scratch, "linearAnionPi*.log"))
merge_logs("logs/planarAnionPi.log", join(scratch, "planarAnionPi*.log"))
merge_logs("logs/methylPi.log", join(scratch, "methylPi*.log"))
merge_logs("logs/additionalInfo.log", join(scratch, "additionalInfo*.log"))
merge_progress_files()

time_stop = time.time()
time_file = open("logs/timeStop.log", "w")
time_file.write(str(time_stop))
time_file.close()

#######profiler stop
# pr.disable()
# s = io.StringIO()
# sortby = 'cumulative'
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()
#
# logFile = open("simpleRun.profile", 'w')
# logFile.write(s.getvalue())
# logFile.close()
##############
