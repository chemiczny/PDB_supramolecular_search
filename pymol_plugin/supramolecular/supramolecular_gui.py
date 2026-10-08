"""
Created on Sun Jul 22 19:39:10 2018

@author: michal
"""

import pandas as pd
from os import path
from time import time

try:
    from pymol import cmd
except ImportError:
    pass

import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk

from cgo_arrow import cgo_arrow


def is_str_int(str2check):
    try:
        int(str2check)
        return True
    except (TypeError, ValueError):
        return False


class SupramolecularGUI:
    def __init__(self, page, parallel_select_function, name):
        self.name = name
        self.page = page
        self.log_data = {
            "logFile": False,
            "data": None,
            "cifDir": None,
            "displaying": False,
            "displayingAround": False,
        }
        self.numerical_parameters = {}
        self.checkbox_vars = {}
        self.list_parameters = {}
        self.sorting_keys2header = {}
        self.unique_keys2header = {}
        self.unique_keys_order = []
        self.interaction_buttons = []
        self.headers = []
        self.current_molecule = {"PdbCode": None}
        self.additional_checkboxes = []
        self.data_is_merged = False
        self.arrow_name = "DefaultArrow"
        self.arrow_color = "blue red"

        self.parallel_selection = False
        self.parallel_select_function = parallel_select_function
        self.external_selection = False

        self.num_of_sorting_menu = 2

        self.actual_displaying = {"rowData": None, "selectionTime": -1}

        self.recording = False
        self.file2record = ""

    def grid(self):
        self.grid_log_file_part()
        self.grid_numerical_parameters()
        self.grid_around_part()
        self.grid_additional_checkboxes()
        self.grid_list_parameters()
        self.grid_sorting_parameters()
        self.grid_unique_parameters()
        self.grid_search_widgets()
        self.grid_tree()
        self.grid_save_filtered()
        self.grid_interaction_buttons()

    def grid_log_file_part(self):
        self.but_log = tk.Button(
            self.page, text="Load log file", command=self.get_log_file, width=10
        )
        self.but_log.grid(row=0, column=0)

        self.ent_log = tk.Entry(self.page, width=20)
        self.ent_log.configure(state="readonly")
        self.ent_log.grid(row=0, column=1, columnspan=2)

        self.lab_use = tk.Label(self.page, text="Use filter")
        self.lab_use.grid(row=0, column=3)

    def get_log_file(self):
        self.log_data["logFile"] = filedialog.askopenfilename(
            title="Select file",
            filetypes=(
                ("Log files", "*.log"),
                ("Txt files", "*.txt"),
                ("CSV files", "*.csv"),
                ("Dat files", "*.dat"),
                ("all files", "*.*"),
            ),
        )
        self.open_log_file()

    def open_log_file(self):
        if not self.log_data["logFile"]:
            return

        self.ent_log.configure(state="normal")
        self.ent_log.delete(0, "end")
        self.ent_log.insert(0, self.log_data["logFile"])
        self.ent_log.configure(state="readonly")

        try:
            self.log_data["data"] = pd.read_csv(
                self.log_data["logFile"], sep="\t"
            ).fillna("NA")
            self.update_menu()

            if self.recording:
                generated_script = open(self.file2record, "a")

                generated_script.write(
                    self.name
                    + ' = pd.read_csv( "'
                    + self.log_data["logFile"]
                    + '", sep = "\\t").fillna("NA") \n'
                )

                generated_script.close()

        except Exception:
            messagebox.showwarning(
                title="Error!", message="Pandas cannot parse this file"
            )

    def update_menu(self):
        for parameter in self.list_parameters:
            parameters_data = self.log_data["data"][
                self.list_parameters[parameter]["header"]
            ].unique()
            parameters_data = sorted(parameters_data)
            self.log_data[parameter] = parameters_data

            self.list_parameters[parameter]["listbox"].delete(0, "end")
            for row in parameters_data:
                self.list_parameters[parameter]["listbox"].insert("end", row)

    def set_numerical_parameters(self, numerical_parameters):
        self.numerical_parameters = numerical_parameters

    def grid_numerical_parameters(self):
        self.actual_row = 1

        for parameter in sorted(self.numerical_parameters.keys()):
            self.numerical_parameters[parameter]["entry_low"] = tk.Entry(
                self.page, width=10
            )
            self.numerical_parameters[parameter]["entry_low"].grid(
                row=self.actual_row, column=0
            )

            self.numerical_parameters[parameter]["label"] = tk.Label(
                self.page, width=10, text="< " + parameter + " <"
            )
            self.numerical_parameters[parameter]["label"].grid(
                row=self.actual_row, column=1
            )

            self.numerical_parameters[parameter]["entry_high"] = tk.Entry(
                self.page, width=10
            )
            self.numerical_parameters[parameter]["entry_high"].grid(
                row=self.actual_row, column=2
            )

            self.checkbox_vars[parameter] = tk.IntVar()

            self.numerical_parameters[parameter]["checkbox"] = tk.Checkbutton(
                self.page, variable=self.checkbox_vars[parameter]
            )
            self.numerical_parameters[parameter]["checkbox"].grid(
                row=self.actual_row, column=3
            )

            self.actual_row += 1

    def list_filter(self, key):
        if not self.log_data["logFile"]:
            return

        template = self.list_parameters[key]["entry"].get()
        template = template.upper()

        template_len = len(template)
        self.list_parameters[key]["listbox"].delete(0, "end")

        for row in self.log_data[key]:
            if template == row[:template_len]:
                self.list_parameters[key]["listbox"].insert("end", row)

    def grid_around_part(self):
        self.lab_around = tk.Label(self.page, text="Around")
        self.lab_around.grid(row=self.actual_row, column=0)

        self.chkvar_around = tk.IntVar()

        self.chk_around = tk.Checkbutton(
            self.page, variable=self.chkvar_around, command=self.show_around
        )
        self.chk_around.grid(row=self.actual_row, column=1)

        self.ent_around = tk.Entry(self.page, width=10)
        self.ent_around.grid(row=self.actual_row, column=2)
        self.ent_around.insert("end", "5.0")

        self.actual_row += 1

    def show_around(self):
        if not self.log_data["displaying"]:
            return

        radius = self.ent_around.get()
        try:
            radius = float(radius)
        except ValueError:
            return

        if self.chkvar_around.get() > 0 and not self.log_data["displayingAround"]:
            selection_around_name = "suprAround"
            cmd.select(
                selection_around_name,
                "byres ( suprSelection around " + str(radius) + " ) ",
            )
            cmd.show("lines", selection_around_name)
            self.log_data["displayingAround"] = True
        elif self.chkvar_around.get() == 0 and self.log_data["displayingAround"]:
            cmd.hide("lines", "suprAround")
            cmd.delete("suprAround")
            self.log_data["displayingAround"] = False

        cmd.deselect()

    def set_list_parameters(self, list_parameters):
        self.list_parameters = list_parameters

    def grid_list_parameters(self):
        self.actual_column = 4
        for parameter in sorted(self.list_parameters.keys()):
            self.list_parameters[parameter]["label"] = tk.Label(
                self.page, text=parameter
            )
            self.list_parameters[parameter]["label"].grid(
                row=0, column=self.actual_column
            )

            self.checkbox_vars[parameter] = tk.IntVar()

            self.list_parameters[parameter]["checkbox"] = tk.Checkbutton(
                self.page, variable=self.checkbox_vars[parameter]
            )
            self.list_parameters[parameter]["checkbox"].grid(
                row=0, column=self.actual_column + 1
            )

            self.list_parameters[parameter]["listbox"] = tk.Listbox(
                self.page, width=10, height=8, exportselection=False
            )
            self.list_parameters[parameter]["listbox"].grid(
                row=0, column=self.actual_column, rowspan=8, columnspan=2
            )
            self.list_parameters[parameter]["listbox"].bind(
                "<<ListboxSelect>>",
                lambda e, arg=parameter: self.move_selected_from_listbox(e, arg),
            )

            self.list_parameters[parameter]["entry"] = tk.Entry(self.page, width=5)
            self.list_parameters[parameter]["entry"].grid(
                row=7, column=self.actual_column
            )

            self.list_parameters[parameter]["button"] = tk.Button(
                self.page,
                width=2,
                text="*",
                command=lambda arg=parameter: self.list_filter(arg),
            )
            self.list_parameters[parameter]["button"].grid(
                row=7, column=self.actual_column + 1
            )

            self.list_parameters[parameter]["listboxSelected"] = tk.Listbox(
                self.page, width=10, height=6, exportselection=False
            )
            self.list_parameters[parameter]["listboxSelected"].grid(
                row=8, column=self.actual_column, rowspan=8, columnspan=2
            )

            self.list_parameters[parameter]["buttonSelClear"] = tk.Button(
                self.page,
                width=2,
                text="clear",
                command=lambda arg=parameter: self.clear_listbox_selected(arg),
            )
            self.list_parameters[parameter]["buttonSelClear"].grid(
                row=16, column=self.actual_column
            )

            self.list_parameters[parameter]["buttonSelDel"] = tk.Button(
                self.page,
                width=2,
                text="del",
                command=lambda arg=parameter: self.remove_from_list_box_selected(arg),
            )
            self.list_parameters[parameter]["buttonSelDel"].grid(
                row=16, column=self.actual_column + 1
            )

            self.actual_column += 2

    def grid_interaction_buttons(self):
        if not self.interaction_buttons:
            return

        current_col = 18
        for button in self.interaction_buttons:
            name = button["name"]

            new_button = tk.Button(
                self.page,
                text=name,
                command=lambda arg=button["headers"]: self.show_more_interactions(arg),
            )
            new_button.grid(row=25, column=current_col, columnspan=2)
            current_col += 2

    def move_selected_from_listbox(self, event, key):
        list_index = self.list_parameters[key]["listbox"].curselection()
        if list_index:
            selection = str(self.list_parameters[key]["listbox"].get(list_index))
            already_selected = self.list_parameters[key]["listboxSelected"].get(
                0, "end"
            )

            if selection not in already_selected:
                self.list_parameters[key]["listboxSelected"].insert("end", selection)

    def clear_listbox_selected(self, key):
        self.list_parameters[key]["listboxSelected"].delete(0, "end")

    def remove_from_list_box_selected(self, key):
        list_index = self.list_parameters[key]["listboxSelected"].curselection()
        if list_index:
            self.list_parameters[key]["listboxSelected"].delete(list_index, list_index)

    def set_sorting_parameters(self, keys2header, keys_col):
        self.sorting_keys2header = keys2header
        self.sorting_keys_col1 = keys_col
        self.sorting_keys_col2 = ["Ascd", "Desc"]

        self.sorting_menu = []

    def grid_sorting_parameters(self):
        if not self.sorting_keys2header:
            return

        for i in range(self.num_of_sorting_menu):
            self.sorting_menu.append({})

            self.sorting_menu[i]["label"] = tk.Label(self.page, text="Sorting" + str(i))
            self.sorting_menu[i]["label"].grid(row=0, column=self.actual_column)

            self.sorting_menu[i]["chk_value"] = tk.IntVar()

            self.sorting_menu[i]["chk_butt"] = tk.Checkbutton(
                self.page, variable=self.sorting_menu[i]["chk_value"]
            )
            self.sorting_menu[i]["chk_butt"].grid(row=0, column=self.actual_column + 1)

            self.sorting_menu[i]["sorting_key"] = tk.IntVar()

            actual_row = 1
            for value, key in enumerate(self.sorting_keys_col1):
                self.sorting_menu[i][key] = tk.Radiobutton(
                    self.page,
                    text=key,
                    variable=self.sorting_menu[i]["sorting_key"],
                    value=value,
                    indicatoron=0,
                    width=8,
                )
                self.sorting_menu[i][key].grid(
                    row=actual_row, column=self.actual_column, columnspan=2
                )
                actual_row += 1

            self.sorting_menu[i]["sortingTypeValue"] = tk.IntVar()

            for value, key in enumerate(self.sorting_keys_col2):
                self.sorting_menu[i][key] = tk.Radiobutton(
                    self.page,
                    text=key,
                    variable=self.sorting_menu[i]["sortingTypeValue"],
                    value=value,
                    indicatoron=0,
                    width=8,
                )
                self.sorting_menu[i][key].grid(
                    row=actual_row, column=self.actual_column, columnspan=2
                )
                actual_row += 1

            self.actual_column += 2

    def set_unique_parameters(self, new_unique_parameters, new_unique_keys_order):
        self.unique_keys2header = new_unique_parameters
        self.unique_keys_order = new_unique_keys_order

    def grid_unique_parameters(self):
        if not self.unique_keys2header:
            return

        unique_label = tk.Label(self.page, text="Unique")
        unique_label.grid(row=0, column=self.actual_column, columnspan=2)

        actual_row = 1

        self.unique_menu = []

        for i, key in enumerate(self.unique_keys_order):
            self.unique_menu.append({})

            new_label = tk.Label(self.page, text=key)
            new_label.grid(row=actual_row, column=self.actual_column)

            self.unique_menu[i]["headers"] = self.unique_keys2header[key]
            self.unique_menu[i]["chk_value"] = tk.IntVar()

            new_chk = tk.Checkbutton(
                self.page, variable=self.unique_menu[i]["chk_value"]
            )
            new_chk.grid(row=actual_row, column=self.actual_column + 1)

            actual_row += 1

        count_button = tk.Button(
            self.page, text="Count", width=8, command=self.count_unique
        )
        count_button.grid(row=actual_row, column=self.actual_column, columnspan=2)

        actual_row += 1

        self.unique_count_entry = tk.Entry(self.page, width=8)
        self.unique_count_entry.grid(
            row=actual_row, column=self.actual_column, columnspan=2
        )
        self.unique_count_entry.configure(state="readonly")

        actual_row += 1

        leave_only_button = tk.Button(
            self.page, text="Leave only", width=8, command=self.leave_only_unique
        )
        leave_only_button.grid(row=actual_row, column=self.actual_column, columnspan=2)

        self.actual_column += 2

    def get_selected_unique_headers(self):
        unique_headers = []

        for unique_data in self.unique_menu:
            if unique_data["chk_value"].get() > 0:
                unique_headers += unique_data["headers"]

        return unique_headers

    def count_unique(self):
        unique_headers = self.get_selected_unique_headers()

        if "filtered" not in self.log_data:
            return

        unique_data = self.log_data["filtered"].drop_duplicates(subset=unique_headers)
        row_number = unique_data.shape[0]

        self.unique_count_entry.configure(state="normal")
        self.unique_count_entry.delete(0, "end")
        self.unique_count_entry.insert(0, str(row_number))
        self.unique_count_entry.configure(state="readonly")

    def leave_only_unique(self):
        unique_headers = self.get_selected_unique_headers()

        if "filtered" not in self.log_data:
            return

        unique_data = self.log_data["filtered"].drop_duplicates(subset=unique_headers)
        self.print_filter_results(unique_data)

    def grid_save_filtered(self):
        self.but_save_filtered = tk.Button(
            self.page, width=7, command=self.save_filtered, text="Save filtered"
        )
        self.but_save_filtered.grid(row=25, column=16, columnspan=2)

    def set_additional_checkboxes(self, additional_checkboxes):
        self.additional_checkboxes = additional_checkboxes

    def grid_additional_checkboxes(self):
        if self.additional_checkboxes:
            for i in range(len(self.additional_checkboxes)):
                self.additional_checkboxes[i]["labTk"] = tk.Label(
                    self.page, text=self.additional_checkboxes[i]["label"]
                )
                self.additional_checkboxes[i]["labTk"].grid(
                    row=self.actual_row, column=0
                )

                self.additional_checkboxes[i]["chkVar"] = tk.IntVar()

                self.additional_checkboxes[i]["chk"] = tk.Checkbutton(
                    self.page, variable=self.additional_checkboxes[i]["chkVar"]
                )
                self.additional_checkboxes[i]["chk"].grid(row=self.actual_row, column=1)

                self.actual_row += 1

    def save_filtered(self):
        if "filtered" not in self.log_data:
            return

        file2save = filedialog.asksaveasfilename(
            defaultextension=".log",
            filetypes=(
                ("Log files", "*.log"),
                ("Txt files", "*.txt"),
                ("CSV files", "*.csv"),
                ("Dat files", "*.dat"),
                ("all files", "*.*"),
            ),
        )
        if file2save:
            self.log_data["filtered"].to_csv(file2save, sep="\t")

        if file2save and self.recording:
            generated_script = open(self.file2record, "a")

            generated_script.write(
                self.name + '_temp.to_csv( "' + file2save + '", sep = "\\t")\n'
            )

            generated_script.close()

    def apply_filter(self):
        if self.log_data["logFile"] is False:
            messagebox.showwarning(title="Warning", message="Log file not selected")
            return

        actual_data = self.log_data["data"]

        if self.recording:
            generated_script = open(self.file2record, "a")

            generated_script.write("\n" + self.name + "_temp = " + self.name + "\n")

            generated_script.close()

        for key in self.checkbox_vars:
            if self.checkbox_vars[key].get() > 0:
                if key in self.list_parameters:
                    selected_values = self.list_parameters[key]["listboxSelected"].get(
                        0, "end"
                    )
                    if selected_values:
                        actual_data = actual_data[
                            actual_data[self.list_parameters[key]["header"]]
                            .astype(str)
                            .isin(selected_values)
                        ]

                    if selected_values and self.recording:
                        generated_script = open(self.file2record, "a")

                        temp_name = self.name + "_temp"
                        selected_values_str = str(list(selected_values))
                        generated_script.write(
                            temp_name
                            + " = "
                            + temp_name
                            + "[ "
                            + temp_name
                            + '[ "'
                            + str(self.list_parameters[key]["header"])
                            + '"].astype(str).isin('
                            + selected_values_str
                            + " )]\n"
                        )

                        generated_script.close()

                elif key in self.numerical_parameters:
                    min_value = self.numerical_parameters[key]["entry_low"].get()
                    max_value = self.numerical_parameters[key]["entry_high"].get()

                    try:
                        min_value = float(min_value)
                        actual_data = actual_data[
                            actual_data[self.numerical_parameters[key]["header"]]
                            > min_value
                        ]

                        if self.recording:
                            generated_script = open(self.file2record, "a")
                            temp_name = self.name + "_temp"

                            generated_script.write(
                                temp_name
                                + " = "
                                + temp_name
                                + "[ "
                                + temp_name
                                + ' [ "'
                                + self.numerical_parameters[key]["header"]
                                + '" ] > '
                                + str(min_value)
                                + " ] \n"
                            )

                            generated_script.close()
                    except Exception:
                        pass

                    try:
                        max_value = float(max_value)
                        actual_data = actual_data[
                            actual_data[self.numerical_parameters[key]["header"]]
                            < max_value
                        ]

                        if self.recording:
                            generated_script = open(self.file2record, "a")
                            temp_name = self.name + "_temp"

                            generated_script.write(
                                temp_name
                                + " = "
                                + temp_name
                                + "[ "
                                + temp_name
                                + ' [ "'
                                + self.numerical_parameters[key]["header"]
                                + '" ] < '
                                + str(max_value)
                                + " ] \n"
                            )

                            generated_script.close()

                    except Exception:
                        pass

        for ad_chk in self.additional_checkboxes:
            if ad_chk["chkVar"].get() > 0:
                actual_data = ad_chk["func"](actual_data)

                if self.recording:
                    generated_script = open(self.file2record, "a")
                    temp_name = self.name + "_temp"

                    generated_script.write(
                        temp_name
                        + " = "
                        + ad_chk["func"].__name__
                        + "("
                        + temp_name
                        + ")\n"
                    )

                    generated_script.close()

        self.data_is_merged = False
        self.print_filter_results(actual_data)

        if self.recording:
            generated_script = open(self.file2record, "a")

            generated_script.write("\n")

            generated_script.close()

    def print_filter_results(self, actual_data):
        records_found = str(len(actual_data))

        any_sort = False
        columns = []
        ascending = []
        for sort_data in self.sorting_menu:
            to_sort = sort_data["chk_value"].get()

            if to_sort > 0:
                any_sort = True
                item_ind = sort_data["sorting_key"].get()
                sorting_key = self.sorting_keys_col1[item_ind]
                header = self.sorting_keys2header[sorting_key]
                if header in columns:
                    continue
                columns.append(header)

                ascending_actual = sort_data["sortingTypeValue"].get()
                if ascending_actual == 0:
                    ascending.append(True)
                else:
                    ascending.append(False)

        self.ent_records_found.configure(state="normal")
        self.ent_records_found.delete(0, "end")
        self.ent_records_found.insert(0, str(records_found))
        self.ent_records_found.configure(state="readonly")
        if any_sort:
            actual_data = actual_data.sort_values(by=columns, ascending=ascending)

        if any_sort and self.recording:
            generated_script = open(self.file2record, "a")

            temp_name = self.name + "_temp"
            ascending_str = [str(val) for val in ascending]
            ascending_str = " [ " + " , ".join(ascending_str) + " ]"
            generated_script.write(
                temp_name
                + " = "
                + temp_name
                + ".sort_values( by = "
                + str(columns)
                + " , ascending = "
                + ascending_str
                + ")\n"
            )

            generated_script.close()
        self.log_data["filtered"] = actual_data

        start, stop = self._get_range()

        diff = stop - start
        new_start = 0
        new_stop = diff

        self.ent_range_start.delete(0, "end")
        self.ent_range_start.insert("end", str(new_start))
        self.ent_range_stop.delete(0, "end")
        self.ent_range_stop.insert("end", str(new_stop))

        self.show_range()
        self.actual_displaying = {"rowData": None, "selectionTime": -1}

    def show_more_interactions(self, headers):
        if "filtered" not in self.log_data:
            return
        current_sel = self.tree_data.focus()
        if current_sel == "":
            return

        row_id = self.tree_data.item(current_sel)["values"][0]
        data = self.log_data["filtered"].iloc[[row_id]]
        pdb_code = data["PDB Code"].values[0]

        if (
            self.current_molecule["PdbCode"]
            and self.current_molecule["PdbCode"] != pdb_code
        ):
            cmd.delete(self.current_molecule["PdbCode"])

        actual_objects = cmd.get_object_list()

        for obj in actual_objects:
            if obj.upper() != pdb_code.upper():
                cmd.delete(obj)

        if self.current_molecule["PdbCode"] != pdb_code:
            if self.log_data["cifDir"] is not None:
                potential_paths = [
                    path.join(self.log_data["cifDir"], pdb_code.lower() + ".cif"),
                    path.join(self.log_data["cifDir"], pdb_code.upper() + ".cif"),
                ]
                cif_found = False
                for file_path in potential_paths:
                    if path.isfile(file_path):
                        cmd.load(file_path)
                        cif_found = True
                        break
                if not cif_found:
                    cmd.fetch(pdb_code)
            else:
                cmd.fetch(pdb_code)

        frame = int(data["Model No"].values[0])
        if frame != 0:
            cmd.frame(frame + 1)

        cmd.hide("everything")
        interactions2print = self.log_data["filtered"]
        for header in headers:
            interactions2print = interactions2print[
                interactions2print[header] == data[header].values[0]
            ]

        selection = ""
        self.delete_arrows()

        for index, row in interactions2print.iterrows():
            arrow_begin, arrow_end = self.get_arrow_from_row(row)
            unique_arrow_name = self.arrow_name + "A" + str(index)
            cgo_arrow(
                arrow_begin,
                arrow_end,
                0.1,
                name=unique_arrow_name,
                color=self.arrow_color,
            )

            if selection == "":
                selection = "(" + self.get_selection_from_row(row) + " ) "
            else:
                selection += "or (" + self.get_selection_from_row(row) + " ) "

        selection_name = "suprSelection"
        cmd.select(selection_name, selection)
        cmd.show("sticks", selection_name)
        cmd.center(selection_name)
        cmd.zoom(selection_name)

        if self.chkvar_around.get() > 0:
            selection_around_name = "suprAround"
            radius = self.ent_around.get()
            try:
                radius = float(radius)
            except ValueError:
                return
            cmd.select(
                selection_around_name,
                "byres ( suprSelection around " + str(radius) + " ) ",
            )
            cmd.show("lines", selection_around_name)
            self.log_data["displayingAround"] = True
        else:
            self.log_data["displayingAround"] = False

        self.log_data["displaying"] = True
        self.current_molecule["PdbCode"] = pdb_code

        cmd.deselect()

        self.actual_displaying["rowData"] = data
        self.actual_displaying["selectionTime"] = time()

    def show_interactions(self):
        if "filtered" not in self.log_data:
            return
        current_sel = self.tree_data.focus()
        if current_sel == "":
            return

        row_id = self.tree_data.item(current_sel)["values"][0]
        data = self.log_data["filtered"].iloc[[row_id]]
        pdb_code = data["PDB Code"].values[0]

        if (
            self.current_molecule["PdbCode"]
            and self.current_molecule["PdbCode"] != pdb_code
        ):
            cmd.delete(self.current_molecule["PdbCode"])

        actual_objects = cmd.get_object_list()

        for obj in actual_objects:
            if obj.upper() != pdb_code.upper():
                cmd.delete(obj)

        if self.current_molecule["PdbCode"] != pdb_code:
            if self.log_data["cifDir"] is not None:
                potential_paths = [
                    path.join(self.log_data["cifDir"], pdb_code.lower() + ".cif"),
                    path.join(self.log_data["cifDir"], pdb_code.upper() + ".cif"),
                ]
                cif_found = False
                for file_path in potential_paths:
                    if path.isfile(file_path):
                        cmd.load(file_path)
                        cif_found = True
                        break
                if not cif_found:
                    cmd.fetch(pdb_code)
            else:
                cmd.fetch(pdb_code)

        frame = int(data["Model No"].values[0])
        if frame != 0:
            cmd.frame(frame + 1)

        selection = self.get_selection(data)

        selection_name = "suprSelection"

        cmd.select(selection_name, selection)
        cmd.show("sticks", selection_name)
        cmd.center(selection_name)
        cmd.zoom(selection_name)

        if self.chkvar_around.get() > 0:
            selection_around_name = "suprAround"
            radius = self.ent_around.get()
            try:
                radius = float(radius)
            except ValueError:
                return
            cmd.select(
                selection_around_name,
                "byres ( suprSelection around " + str(radius) + " ) ",
            )
            cmd.select(selection_around_name, "byres ( suprSelection around 5 ) ")
            cmd.show("lines", selection_around_name)
            self.log_data["displayingAround"] = True
        else:
            self.log_data["displayingAround"] = False

        self.delete_arrows()

        arrow_begin, arrow_end = self.get_arrow(data)
        cgo_arrow(
            arrow_begin, arrow_end, 0.1, name=self.arrow_name, color=self.arrow_color
        )

        self.log_data["displaying"] = True
        self.current_molecule["PdbCode"] = pdb_code

        cmd.deselect()

        self.actual_displaying["rowData"] = data
        self.actual_displaying["selectionTime"] = time()

    def delete_arrows(self):
        for arrow in cmd.get_names_of_type("object:cgo"):
            if "rrow" in arrow:
                cmd.delete(arrow)

    def set_selection_func(self, selection_func):
        self.get_selection = selection_func

    def set_arrow_func(self, arrow_func):
        self.get_arrow = arrow_func

    def show_range(self):
        start, stop = self._get_range()
        if start == stop:
            return

        self._show_range(start, stop)

    def _get_range(self):
        if "filtered" not in self.log_data:
            return 0, 0

        start = self.ent_range_start.get()
        stop = self.ent_range_stop.get()

        if not is_str_int(start) or not is_str_int(stop):
            return 0, 1000

        start = int(start)
        stop = int(stop)

        if stop < start:
            return 0, 1000

        if start < 0:
            return 0, 1000

        return start, stop

    def set_row2values(self, row2values):
        self.get_values = row2values

    def _show_range(self, start, stop):
        self.tree_data.delete(*self.tree_data.get_children())

        row_id = 0
        actual_data = self.log_data["filtered"]
        for index, row in actual_data.iterrows():
            if row_id >= start and row_id < stop:
                self.tree_data.insert("", "end", values=self.get_values(row_id, row))
            row_id += 1
            if row_id >= stop:
                break

    def show_next(self):
        start, stop = self._get_range()
        if start == stop:
            return

        diff = stop - start
        new_start = stop
        new_stop = stop + diff

        data_len = self.log_data["filtered"].shape[0]
        if new_stop > data_len:
            new_stop = data_len
            new_start = data_len - diff

        self.ent_range_start.delete(0, "end")
        self.ent_range_start.insert("end", str(new_start))
        self.ent_range_stop.delete(0, "end")
        self.ent_range_stop.insert("end", str(new_stop))

        self._show_range(new_start, new_stop)

    def show_prev(self):
        start, stop = self._get_range()
        if start == stop:
            return

        diff = stop - start
        new_start = start - diff
        new_stop = start

        if new_start < 0:
            new_stop = diff
            new_start = 0

        self.ent_range_start.delete(0, "end")
        self.ent_range_start.insert("end", str(new_start))
        self.ent_range_stop.delete(0, "end")
        self.ent_range_stop.insert("end", str(new_stop))
        self._show_range(new_start, new_stop)

    def grid_search_widgets(self):
        self.but_apply = tk.Button(
            self.page, width=10, command=self.apply_filter, text="Search"
        )
        self.but_apply.grid(row=22, column=0)

        self.lab_data = tk.Label(self.page, width=10, text="Records found")
        self.lab_data.grid(row=25, column=0)

        self.ent_records_found = tk.Entry(self.page, width=20)
        self.ent_records_found.configure(state="readonly")
        self.ent_records_found.grid(row=25, column=1, columnspan=2)

        self.lab_range = tk.Label(self.page, width=5, text="Range")
        self.lab_range.grid(row=25, column=3)

        self.ent_range_start = tk.Entry(self.page, width=8)
        self.ent_range_start.grid(row=25, column=4, columnspan=2)
        self.ent_range_start.insert("end", 0)

        self.ent_range_stop = tk.Entry(self.page, width=8)
        self.ent_range_stop.grid(row=25, column=6, columnspan=2)
        self.ent_range_stop.insert("end", 100)

        self.but_show_interaction = tk.Button(
            self.page, width=8, command=self.show_interactions, text="Show interact"
        )
        self.but_show_interaction.grid(row=25, column=14, columnspan=2)

        self.but_range_show = tk.Button(
            self.page, width=6, text="Show", command=self.show_range
        )
        self.but_range_show.grid(row=25, column=8, columnspan=2)

        self.but_range_next = tk.Button(
            self.page, width=6, text="Next", command=self.show_next
        )
        self.but_range_next.grid(row=25, column=10, columnspan=2)

        self.but_range_prev = tk.Button(
            self.page, width=6, text="Prev", command=self.show_prev
        )
        self.but_range_prev.grid(row=25, column=12, columnspan=2)

    def set_tree_data(self, tree_data):
        self.headers = tree_data

    def grid_tree(self):
        if not self.headers:
            return

        self.tree_data = ttk.Treeview(
            self.page, columns=self.headers, show="headings", heigh=15
        )
        self.tree_data.bind("<<TreeviewSelect>>", self.tree_parallel_selection)
        for header in self.headers:
            self.tree_data.heading(header, text=header)
            self.tree_data.column(header, width=70)
        self.tree_data.grid(row=30, column=0, columnspan=40)

    def tree_parallel_selection(self, event):
        if not self.data_is_merged:
            return

        if not self.parallel_selection:
            return

        if self.external_selection:
            self.external_selection = False
            return

        current_sel = event.widget.focus()

        row_id = self.tree_data.item(current_sel)["values"][0]
        data = self.log_data["filtered"].iloc[[row_id]]

        self.parallel_select_function(self.name, data)

    def select_row_in_tree(self, value_dict):
        actual_data = self.log_data["filtered"]

        selected_values = actual_data

        for head in value_dict:
            selected_values = selected_values[selected_values[head] == value_dict[head]]

        row2select = selected_values.iloc[0]
        selected_row_index = actual_data.index.get_loc(row2select.name)

        start, stop = self._get_range()
        if selected_row_index < start or selected_row_index > stop:
            diff = stop - start
            new_start = selected_row_index
            new_stop = new_start + diff

            data_len = self.log_data["filtered"].shape[0]
            if new_stop > data_len:
                new_stop = data_len
                new_start = data_len - diff

            self.ent_range_start.delete(0, "end")
            self.ent_range_start.insert("end", str(new_start))
            self.ent_range_stop.delete(0, "end")
            self.ent_range_stop.insert("end", str(new_stop))

            self._show_range(new_start, new_stop)

            start = new_start
            stop = new_stop

        treeview_len = len(self.tree_data.get_children())
        item_id = self.tree_data.get_children()[selected_row_index - start]
        self.tree_data.selection_set(item_id)
        fraction = float(selected_row_index - start) / treeview_len
        self.tree_data.yview_moveto(fraction)

        self.external_selection = True

    def get_state(self):
        state = {
            "numerical_parameters": {},
            "checkboxes": {},
            "additionalCheckboxes": {},
            "listParameters": {},
        }

        for label in self.numerical_parameters:
            state["numerical_parameters"][label] = {}

            state["numerical_parameters"][label]["entry_low"] = (
                self.numerical_parameters[label]["entry_low"].get()
            )
            state["numerical_parameters"][label]["entry_high"] = (
                self.numerical_parameters[label]["entry_high"].get()
            )

        for label in self.checkbox_vars:
            state["checkboxes"][label] = self.checkbox_vars[label].get()

        for label in self.list_parameters:
            state["listParameters"][label] = self.list_parameters[label][
                "listboxSelected"
            ].get(0, "end")

        for obj in self.additional_checkboxes:
            state["additionalCheckboxes"][obj["label"]] = obj["chkVar"].get()

        return state

    def load_state(self, state):

        for label in state["numerical_parameters"]:
            if label in self.numerical_parameters:
                new_value_low = state["numerical_parameters"][label]["entry_low"]
                new_value_high = state["numerical_parameters"][label]["entry_high"]

                self.numerical_parameters[label]["entry_low"].delete(0, "end")
                self.numerical_parameters[label]["entry_low"].insert(0, new_value_low)

                self.numerical_parameters[label]["entry_high"].delete(0, "end")
                self.numerical_parameters[label]["entry_high"].insert(0, new_value_high)

        for label in state["checkboxes"]:
            if label in self.checkbox_vars:
                new_value = int(state["checkboxes"][label])
                self.checkbox_vars[label].set(new_value)

        for label in state["listParameters"]:
            if label in self.list_parameters:
                new_values = state["listParameters"][label]

                self.list_parameters[label]["listboxSelected"].delete(0, "end")
                for val in new_values:
                    self.list_parameters[label]["listboxSelected"].insert("end", val)

        for obj in self.additional_checkboxes:
            label = obj["label"]
            if label in state["additionalCheckboxes"]:
                new_value = int(state["additionalCheckboxes"][label])
                obj["chkVar"].set(new_value)

    def start_recording(self, file2record):
        self.file2record = file2record
        self.recording = True

    def stop_recording(self):
        self.recording = False
