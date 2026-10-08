"""
Created on Fri Aug  3 10:31:51 2018

@author: michal
"""

from anion_pi import AnionPiGUI
from pi_pi import PiPiGUI
from cation_pi import CationPiGUI
from anion_cation import AnionCationGUI
from h_bonds import HBondsGUI
from metal_ligand import MetalLigandGUI

from linear_anion_pi import LinearAnionPiGUI
from planar_anion_pi import PlanarAnionPiGUI

from qm_gui import QMGUI
from job_status_gui import JobStatusGUI
from os import path
import pandas as pd
from cgo_arrow import cgo_arrow
import json

try:
    from pymol import cmd
except ImportError:
    pass

import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog


class SupramolecularComposition:
    def __init__(
        self,
        page_anion_pi,
        page_pi_pi,
        page_cation_pi,
        page_anion_cation,
        page_h_bonds,
        page_metal_ligand,
        page_linear_anion_pi,
        page_planar_anion_pi,
        page_qm,
        page_job_status,
    ):
        self.gui_anion_pi = AnionPiGUI(page_anion_pi, self.parallel_select, "AnionPi")
        self.gui_pi_pi = PiPiGUI(page_pi_pi, self.parallel_select, "PiPi")
        self.gui_cation_pi = CationPiGUI(
            page_cation_pi, self.parallel_select, "CationPi"
        )
        self.gui_anion_cation = AnionCationGUI(
            page_anion_cation, self.parallel_select, "AnionCation"
        )
        self.gui_h_bonds = HBondsGUI(page_h_bonds, self.parallel_select, "HBonds")
        self.gui_metal_ligand = MetalLigandGUI(
            page_metal_ligand, self.parallel_select, "MetalLigand"
        )
        self.gui_linear_anion_pi = LinearAnionPiGUI(
            page_linear_anion_pi, self.parallel_select, "LinearAnionPi"
        )
        self.gui_planar_anion_pi = PlanarAnionPiGUI(
            page_planar_anion_pi, self.parallel_select, "PlanarAnionPi"
        )
        self.gui_qm = QMGUI(page_qm)
        self.gui_job_status = JobStatusGUI(page_job_status)

        self.guis = [
            self.gui_anion_pi,
            self.gui_pi_pi,
            self.gui_cation_pi,
            self.gui_anion_cation,
            self.gui_h_bonds,
            self.gui_metal_ligand,
            self.gui_linear_anion_pi,
            self.gui_planar_anion_pi,
        ]

        self.action_labels = [
            "AnionPi",
            "PiPi",
            "CationPi",
            "AnionCation",
            "HBonds",
            "MetalLigand",
            "LinearAnionPi",
            "PlanarAnionPi",
        ]
        self.action_labels2objects = {
            "AnionPi": self.gui_anion_pi,
            "PiPi": self.gui_pi_pi,
            "CationPi": self.gui_cation_pi,
            "AnionCation": self.gui_anion_cation,
            "HBonds": self.gui_h_bonds,
            "MetalLigand": self.gui_metal_ligand,
            "LinearAnionPi": self.gui_linear_anion_pi,
            "PlanarAnionPi": self.gui_planar_anion_pi,
        }

        self.actual_keys = []

        self.recording = False
        self.file2record = ""

        self.anions_as_individuals = tk.IntVar()
        self.residue_rings_as_individuals = tk.IntVar()

    def merge(self, action_menu):
        selected_data = []
        data2use = {}

        for label in self.action_labels:
            use = action_menu[label]["checkValue"].get() > 0
            exclude = action_menu[label]["checkValueExclude"].get() > 0
            if use or exclude:
                if "filtered" not in self.action_labels2objects[label].log_data:
                    self.action_labels2objects[label].data_is_merged = False
                    continue
                else:
                    if use and exclude:
                        self.action_labels2objects[label].data_is_merged = False
                        continue

                    data2use[label] = use
                    selected_data.append(label)
            else:
                self.action_labels2objects[label].data_is_merged = False

        if len(selected_data) < 2:
            return

        unique_data = []
        data_excluded = []
        self.actual_keys = []
        excluded_keys = []

        merging_variable_names = []
        merging_headers = []

        exclude_variable_names = []
        exclude_headers = []

        for key in selected_data:
            headers = self.get_headers_for_gui(key)

            if data2use[key]:
                if len(unique_data) == 0:
                    unique_data = (
                        self.action_labels2objects[key]
                        .log_data["filtered"][headers]
                        .drop_duplicates()
                    )
                elif len(self.action_labels2objects[key].log_data["filtered"]) > 0:
                    unique_data = pd.merge(
                        unique_data,
                        self.action_labels2objects[key].log_data["filtered"][headers],
                        on=list(set(self.actual_keys) & set(headers)),
                    )
                    unique_data = unique_data.drop_duplicates()
                self.actual_keys = list(set(self.actual_keys + headers))

                merging_variable_names.append(
                    self.action_labels2objects[key].name + "_temp"
                )
                merging_headers.append(headers)

            else:
                if len(data_excluded) == 0:
                    data_excluded = (
                        self.action_labels2objects[key]
                        .log_data["filtered"][headers]
                        .drop_duplicates()
                    )
                elif len(self.action_labels2objects[key].log_data["filtered"]) > 0:
                    data_excluded = pd.merge(
                        data_excluded,
                        self.action_labels2objects[key].log_data["filtered"][headers],
                        on=list(set(self.actual_keys) & set(headers)),
                    )
                    data_excluded = data_excluded.drop_duplicates()
                excluded_keys = list(set(excluded_keys + headers))

                exclude_variable_names.append(
                    self.action_labels2objects[key].name + "_temp"
                )
                exclude_headers.append(headers)

        if self.recording:
            generated_script = open(self.file2record, "a")

            generated_script.write(
                "data_frames_to_merge = [ "
                + " , ".join(merging_variable_names)
                + " ]\n"
            )
            generated_script.write(
                "data_frame_merge_headers = " + str(merging_headers) + "\n"
            )

            generated_script.write(
                "data_frames_to_exclude = [ "
                + " , ".join(exclude_variable_names)
                + " ]\n"
            )
            generated_script.write(
                "data_frame_exclude_headers = " + str(exclude_headers) + "\n\n"
            )

            all_variables = merging_variable_names + exclude_variable_names

            generated_script.write(
                "[ "
                + " , ".join(all_variables)
                + " ] = "
                + "simple_merge(data_frames_to_merge, data_frame_merge_headers, "
                + "data_frames_to_exclude, data_frame_exclude_headers)\n"
            )

            generated_script.close()

        if len(data_excluded) > 0:
            merging_keys = list(set(self.actual_keys) & set(excluded_keys))
            sub_merged = pd.merge(
                unique_data, data_excluded, on=merging_keys, how="left", indicator=True
            )
            unique_data = sub_merged[sub_merged["_merge"] == "left_only"]

        for key in selected_data:
            headers = self.get_headers_for_gui(key)

            if len(unique_data) == 0:
                messagebox.showwarning(
                    title="Merging error!", message="No data left after merge!"
                )
                break

            merging_keys = list(set(self.actual_keys) & set(headers))
            temp_data_frame = unique_data[merging_keys].drop_duplicates()
            merged_data = pd.merge(
                self.action_labels2objects[key].log_data["filtered"],
                temp_data_frame,
                on=merging_keys,
            )
            self.action_labels2objects[key].log_data["filtered"] = merged_data
            self.action_labels2objects[key].print_filter_results(merged_data)
            self.action_labels2objects[key].data_is_merged = True

    def show_all(self, show_menu):
        if not self.actual_keys:
            return

        selected_data = []
        last_selection_menu = ""
        last_selection_time = -1

        for label in self.action_labels:
            new_time = self.action_labels2objects[label].actual_displaying[
                "selectionTime"
            ]
            if new_time > last_selection_time:
                last_selection_time = new_time
                last_selection_menu = label

            if show_menu[label]["checkValue"].get() > 0:
                if "filtered" not in self.action_labels2objects[label].log_data:
                    continue
                else:
                    selected_data.append(label)

        if len(selected_data) < 1:
            return

        if last_selection_time < 0:
            return

        if last_selection_menu not in selected_data:
            return

        self.delete_merged_arrows()
        unique_data = self.action_labels2objects[last_selection_menu].actual_displaying[
            "rowData"
        ]

        last_data_headers = self.get_headers_for_gui(last_selection_menu)

        selection = ""
        key_sum = set(last_data_headers)
        for key in selected_data:
            headers = self.get_headers_for_gui(key)

            merging_keys = list(key_sum & set(headers))
            unique_data = pd.merge(
                self.action_labels2objects[key].log_data["filtered"],
                unique_data,
                on=merging_keys,
            )
            key_sum |= set(headers)

        for key in selected_data:
            headers = self.get_headers_for_gui(key)
            merging_keys = list(set(headers) & key_sum)
            temp_data_frame = unique_data[merging_keys].drop_duplicates()
            merged_data = pd.merge(
                self.action_labels2objects[key].log_data["filtered"],
                temp_data_frame,
                on=merging_keys,
            )
            for index, row in merged_data.iterrows():
                arrow_begin, arrow_end = self.action_labels2objects[
                    key
                ].get_arrow_from_row(row)
                unique_arrow_name = (
                    self.action_labels2objects[key].arrow_name + "A" + str(index)
                )
                cgo_arrow(
                    arrow_begin,
                    arrow_end,
                    0.1,
                    name=unique_arrow_name,
                    color=self.action_labels2objects[key].arrow_color,
                )

                if selection == "":
                    selection = (
                        "("
                        + self.action_labels2objects[key].get_selection_from_row(row)
                        + " ) "
                    )
                else:
                    selection += (
                        "or ("
                        + self.action_labels2objects[key].get_selection_from_row(row)
                        + " ) "
                    )

        selection_name = "suprSelectionFull"
        cmd.select(selection_name, selection)
        cmd.show("sticks", selection_name)
        cmd.center(selection_name)
        cmd.zoom(selection_name)

    def delete_merged_arrows(self):
        for arrow in cmd.get_names_of_type("object:cgo"):
            if "rrowA" in arrow:
                cmd.delete(arrow)

    def read_all_logs_from_dir(self, log_dir):
        for gui_key in self.action_labels2objects:
            basename = gui_key[0].lower() + gui_key[1:] + ".log"
            log_file_name = path.join(log_dir, basename)
            if path.isfile(log_file_name):
                self.action_labels2objects[gui_key].log_data["logFile"] = log_file_name
                self.action_labels2objects[gui_key].open_log_file()

    def select_cif_dir(self, cif_dir):
        for gui in self.guis:
            gui.log_data["cifDir"] = cif_dir

    def grid(self):
        for gui in self.guis:
            gui.grid()

        self.gui_qm.grid()
        self.gui_job_status.grid()

    def load_state(self):
        json_file = filedialog.askopenfilename(
            title="Select file",
            filetypes=(("Json files", "*.json"), ("all files", "*.*")),
        )

        if not json_file:
            return

        with open(json_file, "r") as fp:
            state = json.load(fp)

        for label in self.action_labels2objects:
            if label in state:
                self.action_labels2objects[label].load_state(state[label])

        if "QMGUI" in state:
            self.gui_qm.load_state(state["QMGUI"])

        if "JobStatus" in state:
            self.gui_job_status.load_state(state["JobStatus"])

    def save_state(self):
        state = {}
        for label in self.action_labels2objects:
            state[label] = self.action_labels2objects[label].get_state()

        state["QMGUI"] = self.gui_qm.get_state()
        state["JobStatus"] = self.gui_job_status.get_state()

        file2save = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=(("Json files", "*.json"), ("all files", "*.*")),
        )
        if file2save:
            with open(file2save, "w") as fp:
                json.dump(state, fp)

    def set_parallel_selection(self, state):
        if state == 1:
            for gui in self.guis:
                gui.parallel_selection = True
        else:
            for gui in self.guis:
                gui.parallel_selection = False

    def get_headers_for_gui(self, gui_key):
        last_data_headers = ["PDB Code", "Model No"]

        headers_id = {
            "Pi": ["Pi acid Code", "Pi acid chain", "Piacid id"],
            "Anion": ["Anion code", "Anion chain", "Anion id"],
            "Cation": ["Cation code", "Cation chain", "Cation id"],
            "HBonds": ["Anion code", "Anion chain", "Anion id"],
            "Metal": ["Cation code", "Cation chain", "Cation id"],
            "Ligand": ["Anion code", "Anion chain", "Anion id"],
        }

        if self.anions_as_individuals.get() > 0:
            headers_id["Anion"] = [
                "Anion code",
                "Anion chain",
                "Anion id",
                "Anion group id",
            ]
            headers_id["Ligand"] = [
                "Anion code",
                "Anion chain",
                "Anion id",
                "Anion group id",
            ]
            headers_id["HBonds"] = [
                "Anion code",
                "Anion chain",
                "Anion id",
                "Anion group id",
            ]

        if self.residue_rings_as_individuals.get() > 0:
            headers_id["Pi"] = [
                "Pi acid Code",
                "Pi acid chain",
                "Piacid id",
                "CentroidId",
            ]

        for header_key in headers_id:
            if header_key in gui_key:
                last_data_headers += headers_id[header_key]

        if "Anion" in gui_key and "Cation" in gui_key:
            last_data_headers += headers_id["Pi"]

        return last_data_headers

    def parallel_select(self, name, selected_row):
        selected_headers = self.get_headers_for_gui(name)
        for gui_name in self.action_labels2objects:
            if gui_name == name:
                continue

            if not self.action_labels2objects[gui_name].data_is_merged:
                continue

            headers = self.get_headers_for_gui(gui_name)
            common_headers = set(selected_headers) & set(headers)
            header2value = {}

            for head in common_headers:
                header2value[head] = selected_row[head].values[0]

            self.action_labels2objects[gui_name].select_row_in_tree(header2value)

    def start_recording(self):
        file2record = filedialog.asksaveasfilename(
            defaultextension=".py",
            filetypes=(("Python files", "*.py"), ("all files", "*.*")),
        )

        if not file2record:
            return

        self.file2record = file2record
        self.recording = True

        generated_script = open(self.file2record, "w")
        generated_script.write("#Script generated by SupremolecularAnalyser\n")
        generated_script.write("import pandas as pd\n")
        generated_script.write("from simple_filters import *\n\n")
        generated_script.close()

        for gui in self.guis:
            gui.start_recording(self.file2record)

    def stop_recording(self):
        self.recording = False

        for gui in self.guis:
            gui.stop_recording()
