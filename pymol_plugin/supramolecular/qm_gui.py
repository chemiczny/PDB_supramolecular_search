"""
Created on Mon Feb 25 16:21:16 2019

@author: michal
"""

import os
from copy import deepcopy

try:
    from pymol import cmd
except ImportError:
    pass

import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox


def get_all_selection_names():
    return cmd.get_names("selections", 0)


class QMGUI:
    def __init__(self, page):
        self.page = page
        self.keyword_set = {}

        self.charge_total = None
        self.spin_total = None
        self.charge_guest = None
        self.spin_guest = None
        self.charge_host = None
        self.spin_host = None

    def grid_basic_labels(self):
        actual_row = 1
        charge_label = tk.Label(self.page, text="charge")
        charge_label.grid(row=actual_row, column=0)

        actual_row += 1

        spin_label = tk.Label(self.page, text="spin")
        spin_label.grid(row=actual_row, column=0)

        actual_row += 1

        sele_label = tk.Label(self.page, text="sele")
        sele_label.grid(row=actual_row, column=0, rowspan=2)

        actual_row += 2

        frozen_label = tk.Label(self.page, text="frozen")
        frozen_label.grid(row=actual_row, column=0, rowspan=2)

        actual_row += 2

        add_atoms_label = tk.Label(self.page, text="Add atoms")
        add_atoms_label.grid(row=actual_row, column=0, rowspan=2)

        actual_row += 2

        self.count_atoms_button = tk.Button(
            self.page, width=5, text="count at", command=self.count_atoms
        )
        self.count_atoms_button.grid(row=actual_row, column=0)

    def count_atoms(self):
        try:
            model_host = cmd.get_model("host")
            host_atom_no = len(model_host.atom)
        except Exception:
            host_atom_no = 0

        try:
            model_guest = cmd.get_model("guest")
            guest_atom_no = len(model_guest.atom)
        except Exception:
            guest_atom_no = 0

        complex_atom_no = host_atom_no + guest_atom_no

        self.guest_atoms_no.delete(0, "end")
        self.guest_atoms_no.insert("end", str(guest_atom_no))

        self.host_atoms_no.delete(0, "end")
        self.host_atoms_no.insert("end", str(host_atom_no))

        self.complex_atoms_no.delete(0, "end")
        self.complex_atoms_no.insert("end", str(complex_atom_no))

    def grid_host(self):
        actual_row = 0

        host_label = tk.Label(self.page, text="Host")
        host_label.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_charge = tk.Entry(self.page, width=5)
        self.host_charge.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_spin = tk.Entry(self.page, width=5)
        self.host_spin.grid(row=actual_row, column=1)
        self.host_spin.insert("end", "1")

        actual_row += 1

        self.host_sele = tk.Button(
            self.page, width=5, text="get", command=self.read_host
        )
        self.host_sele.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_add_sele = tk.Button(
            self.page, width=5, text="add", command=self.add_sele2_host
        )
        self.host_add_sele.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_frozen_clear = tk.Button(
            self.page, width=5, text="clear", command=self.clear_host_frozen
        )
        self.host_frozen_clear.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_frozen_default = tk.Button(
            self.page, width=5, text="default", command=self.default_host_frozen
        )
        self.host_frozen_default.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_add_h = tk.Button(
            self.page, width=5, text="addH", command=self.add_h_host
        )
        self.host_add_h.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_add_nca = tk.Button(
            self.page, width=5, text="chop", command=self.host_chop
        )
        self.host_add_nca.grid(row=actual_row, column=1)

        actual_row += 1

        self.host_atoms_no = tk.Entry(self.page, width=5)
        self.host_atoms_no.grid(row=actual_row, column=1)

        actual_row += 1

    def add_sele2_host(self):
        try:
            state_no = cmd.get_state()
            cmd.create("host", "%host or sele", state_no)
        except Exception:
            state_no = cmd.get_state()
            cmd.create("host", "sele", state_no)

    def clear_host_frozen(self):
        cmd.select("hostFrozen", "none")

    def default_host_frozen(self):
        cmd.select("hostFrozen", " %hostFrozen or ( host and ( name CA ) ) ")

    def add_h_host(self):
        cmd.h_add("host")

    def host_chop(self):
        self.chop("host")

    def read_host(self):
        try:
            state_no = cmd.get_state()
            cmd.create("host", "sele", state_no)
            cmd.select("hostFrozen", "none")
        except Exception:
            print("lo kurla")

    def read_host_from_sele_around(self):
        try:
            state_no = cmd.get_state()
            radius = self.sele_radius.get()
            cmd.create("host", "byres ( sele around " + radius + ")", state_no)
        except Exception:
            print("lo kurla")

    def grid_guest(self):
        actual_row = 0
        guest_label = tk.Label(self.page, text="Guest")
        guest_label.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_charge = tk.Entry(self.page, width=5)
        self.guest_charge.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_spin = tk.Entry(self.page, width=5)
        self.guest_spin.grid(row=actual_row, column=2)
        self.guest_spin.insert("end", "1")

        actual_row += 1

        self.guest_sele = tk.Button(
            self.page, text="get", width=5, command=self.read_guest
        )
        self.guest_sele.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_sele_add = tk.Button(
            self.page, text="add", width=5, command=self.add_sele2_guest
        )
        self.guest_sele_add.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_frozen_clear = tk.Button(
            self.page, width=5, text="clear", command=self.clear_guest_frozen
        )
        self.guest_frozen_clear.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_frozen_default = tk.Button(
            self.page, width=5, text="default", command=self.default_guest_frozen
        )
        self.guest_frozen_default.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_add_h = tk.Button(
            self.page, width=5, text="addH", command=self.add_h_guest
        )
        self.guest_add_h.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_add_nca = tk.Button(
            self.page, width=5, text="chop", command=self.guest_chop
        )
        self.guest_add_nca.grid(row=actual_row, column=2)

        actual_row += 1

        self.guest_atoms_no = tk.Entry(self.page, width=5)
        self.guest_atoms_no.grid(row=actual_row, column=2)

        actual_row += 1

    def add_sele2_guest(self):
        try:
            state_no = cmd.get_state()
            cmd.create("guest", "%guest or sele", state_no)
        except Exception:
            state_no = cmd.get_state()
            cmd.create("guest", "sele", state_no)

    def chop(self, sele):
        try:
            model = cmd.get_model(sele)
        except Exception:
            return

        chain2res_id2atom_names = {}

        for atom in model.atom:
            resnum = int(atom.resi)
            pdbname = atom.name
            chain = atom.chain

            if chain not in chain2res_id2atom_names:
                chain2res_id2atom_names[chain] = {resnum: set([pdbname])}
            else:
                if resnum in chain2res_id2atom_names[chain]:
                    chain2res_id2atom_names[chain][resnum].add(pdbname)
                else:
                    chain2res_id2atom_names[chain][resnum] = set([pdbname])

        for chain in chain2res_id2atom_names:
            res_id2atom_names = chain2res_id2atom_names[chain]

            res_ids2add_n = set([])
            res_ids2add_c = set([])

            for resnum in res_id2atom_names:
                if res_id2atom_names[resnum] != set(
                    ["CA", "C", "O"]
                ) and res_id2atom_names[resnum] != set(["CA", "N"]):
                    res_ids2add_c.add(resnum + 1)
                    res_ids2add_n.add(resnum - 1)

            res_ids2add_n -= set(res_id2atom_names.keys())
            res_ids2add_c -= set(res_id2atom_names.keys())

            state_no = cmd.get_state()
            for resnum in res_ids2add_c:
                cmd.create(
                    sele,
                    " %"
                    + sele
                    + " or ( ( resi "
                    + str(resnum)
                    + " and chain "
                    + chain
                    + " ) and ( name CA or name N) ) ",
                    state_no,
                )
                cmd.bond(
                    "%" + sele + " and name N and resi " + str(resnum),
                    "%" + sele + " and name C and resi " + str(resnum - 1),
                    1,
                )

            for resnum in res_ids2add_n:
                cmd.create(
                    sele,
                    " %"
                    + sele
                    + " or ( ( resi "
                    + str(resnum)
                    + " and chain "
                    + chain
                    + " ) and ( name CA or name C or name O) ) ",
                    state_no,
                )
                cmd.bond(
                    "%" + sele + " and name C and resi " + str(resnum),
                    "%" + sele + " and name N and resi " + str(resnum + 1),
                    1,
                )

        cmd.show("sticks", sele)

    def guest_chop(self):
        self.chop("guest")

    def clear_guest_frozen(self):
        cmd.select("guestFrozen", "none")

    def default_guest_frozen(self):
        cmd.select("guestFrozen", " %guestFrozen or ( guest and ( name CA  ) ) ")

    def add_h_guest(self):
        cmd.h_add("guest")

    def read_guest(self):
        try:
            state_no = cmd.get_state()
            cmd.create("guest", "sele", state_no)
            cmd.select("guestFrozen", "none")
        except Exception:
            print("lo kurla")

    def read_guest_from_sele_around(self):
        try:
            state_no = cmd.get_state()
            radius = self.sele_radius.get()
            cmd.create("guest", "byres ( sele around " + radius + ")", state_no)
        except Exception:
            print("lo kurla")

    def grid_complex(self):
        complex_label = tk.Label(self.page, text="Complex")
        complex_label.grid(row=0, column=3)

        self.complex_charge = tk.Entry(self.page, width=5)
        self.complex_charge.grid(row=1, column=3)

        self.complex_spin = tk.Entry(self.page, width=5)
        self.complex_spin.grid(row=2, column=3)
        self.complex_spin.insert("end", "1")

        self.complex_atoms_no = tk.Entry(self.page, width=5)
        self.complex_atoms_no.grid(row=9, column=3)

    def grid_write(self):
        write_button = tk.Button(self.page, text="write", command=self.write)
        write_button.grid(row=10, column=0)

    def write(self):
        self.charge_total = self.complex_charge.get()
        self.spin_total = self.complex_spin.get()
        self.charge_guest = self.guest_charge.get()
        self.spin_guest = self.guest_spin.get()
        self.charge_host = self.host_charge.get()
        self.spin_host = self.host_spin.get()

        if not self.charge_total or not self.charge_guest or not self.charge_host:
            messagebox.showwarning(
                title="Error!", message="Please fill the charge data"
            )
            return

        if not self.spin_total or not self.spin_guest or not self.spin_host:
            messagebox.showwarning(title="Error!", message="Please fill the spin data")
            return

        options = {"mustexist": False, "title": "Job directory selection"}
        directory = filedialog.askdirectory(**options)
        if not directory:
            return

        if not os.path.isdir(directory):
            os.makedirs(directory)

        basename = os.path.basename(directory)

        cmd.save(os.path.join(directory, basename + ".pdb"), "host or guest")

        additional = open(os.path.join(directory, basename + ".dat"), "w")
        objects = cmd.get_object_list()
        for obj in objects:
            if obj != "host" and obj != "guest":
                additional.write(obj + "\n")

        additional.close()

        if self.keyword_set:
            for set_name in self.keyword_set:
                self.combine_set(directory, set_name)

        else:
            self.write_input_job(directory, basename, {})

    def combine_set(self, directory, set_name):
        queue = [{}]

        for key in self.keyword_set[set_name]:
            new_queue = []

            for value in self.keyword_set[set_name][key]:
                for element in queue:
                    new_element = deepcopy(element)
                    new_element[key] = value
                    new_queue.append(new_element)

            queue = new_queue

        for job_dict in queue:
            basename = set_name
            for key in job_dict:
                new_part = job_dict[key]
                new_part = new_part.replace(" ", "")
                new_part = new_part.replace("=", "")
                new_part = new_part.replace("(", "")
                new_part = new_part.replace(")", "")

                basename += "_" + new_part

            new_directory = os.path.join(directory, basename)
            os.makedirs(new_directory)

            self.write_input_job(new_directory, basename, job_dict)

    def write_input_job(self, directory, basename, keyword_dict):

        slurm_file = os.path.join(directory, basename + ".slurm")
        inp_file = os.path.join(directory, basename + ".inp")
        inp_file_basename = os.path.basename(inp_file)

        slurm_f = open(slurm_file, "w")

        slurm_head = self.slurm_text_g16.get("1.0", "end")
        slurm_head = self.update_text(slurm_head, keyword_dict)

        slurm_f.write(slurm_head)
        slurm_f.write("\n")

        slurm_f.write("\nmodule add plgrid/apps/gaussian/g16.A.03\n\n")
        slurm_f.write("g16 " + inp_file_basename + "\n\n")

        slurm_f.close()

        inp_f = open(inp_file, "w")

        inp_f.write("%Chk=" + inp_file_basename.replace(".inp", ".chk") + "\n")
        route_section = self.route_section.get("1.0", "end")
        route_section = self.update_text(route_section, keyword_dict)
        inp_f.write(route_section)

        inp_f.write("\nEmilka jest najpiekniejsza!\n\n")

        inp_f.write(
            self.charge_total
            + ","
            + self.spin_total
            + " "
            + self.charge_guest
            + ","
            + self.spin_guest
            + " "
            + self.charge_host
            + ","
            + self.spin_host
            + "\n"
        )

        model = cmd.get_model("guest and guestFrozen")
        self.write_model_to_file(model, inp_f, 1, True)

        model = cmd.get_model("guest and not guestFrozen")
        self.write_model_to_file(model, inp_f, 1)

        model = cmd.get_model("host and hostFrozen")
        self.write_model_to_file(model, inp_f, 2, True)

        model = cmd.get_model("host and not hostFrozen")
        self.write_model_to_file(model, inp_f, 2)

        inp_f.write("\n")

        additional_input = self.additional_section.get("1.0", "end")
        additional_input = self.update_text(additional_input, keyword_dict)
        inp_f.write(additional_input)
        inp_f.write("\n\n")

        inp_f.close()

    def write_model_to_file(self, model, file2append, fragment_no, frozen=False):
        for atom in model.atom:
            resnum = atom.resi
            resname = atom.resn
            chain = atom.chain
            pdbname = atom.name
            element = atom.symbol

            x = str(atom.coord[0])
            y = str(atom.coord[1])
            z = str(atom.coord[2])

            if not frozen:
                file2append.write(
                    element
                    + "(Fragment="
                    + str(fragment_no)
                    + ", PDBName="
                    + pdbname
                    + ", ResName="
                    + resname
                    + ", ResNum="
                    + str(resnum)
                    + "_"
                    + chain
                    + ") "
                    + x
                    + " "
                    + y
                    + " "
                    + z
                    + "\n"
                )
            else:
                file2append.write(
                    element
                    + "(Fragment="
                    + str(fragment_no)
                    + ", PDBName="
                    + pdbname
                    + ", ResName="
                    + resname
                    + ", ResNum="
                    + str(resnum)
                    + "_"
                    + chain
                    + ") -1"
                    + x
                    + " "
                    + y
                    + " "
                    + z
                    + "\n"
                )

    def get_keywords_from_text(self, text):
        keywords = []

        inside_keyword = False

        for letter in text:
            if letter == "{":
                inside_keyword = True
                keywords.append("")
                continue
            elif letter == "}":
                inside_keyword = False

            if inside_keyword:
                keywords[-1] += letter

        return keywords

    def update_text(self, text, key_dict):
        keywords = self.get_keywords_from_text(text)

        for key in keywords:
            if key in key_dict:
                text = text.replace("{" + key + "}", key_dict[key])
            else:
                text = text.replace("{" + key + "}", " ")

        return text

    def grid_gaussian_route_section(self):
        route_section_label = tk.Label(self.page, text="Gaussian route section")
        route_section_label.grid(row=1, column=5, columnspan=5)

        self.route_section = tk.Text(self.page, width=50, height=10)
        self.route_section.grid(row=2, column=5, columnspan=5, rowspan=5)
        self.route_section.insert(
            "end",
            "%Mem=100GB\n"
            "#P B3LYP/6-31G(d,p)\n"
            "# Opt Counterpoise=2\n"
            "# SCRF(Solvent=Water, Read)\n"
            "# Gfinput IOP(6/7=3)  Pop=full  Density  Test \n"
            "# Units(Ang,Deg)",
        )

    def grid_gaussian_additional_input(self):
        additional_input_label = tk.Label(self.page, text="Additional input")
        additional_input_label.grid(row=7, column=5, columnspan=5)

        self.additional_section = tk.Text(self.page, width=50, height=10)
        self.additional_section.grid(row=8, column=5, columnspan=5, rowspan=5)
        self.additional_section.insert("end", "eps=4\n")

    def grid_slurm_section(self):
        slurm_section = tk.Label(self.page, text="Slurm input")
        slurm_section.grid(row=1, column=10, columnspan=10)

        self.slurm_text_g16 = tk.Text(self.page, width=50, height=10)
        self.slurm_text_g16.grid(row=2, column=10, columnspan=5, rowspan=5)
        self.slurm_text_g16.insert(
            "end",
            "#!/bin/env bash\n"
            "#SBATCH --nodes=1\n"
            "#SBATCH --cpus-per-task=24\n"
            "#SBATCH --time=70:00:00\n"
            "##### Nazwa kolejki\n"
            "#SBATCH -p plgrid\n",
        )

    def grid_set_keywords_section(self):
        set_label = tk.Label(self.page, text="Set")
        set_label.grid(row=7, column=10, columnspan=2)

        self.set_listbox = tk.Listbox(
            self.page, width=12, height=9, exportselection=False
        )
        self.set_listbox.grid(row=8, column=10, columnspan=2, rowspan=5)
        self.set_listbox.bind("<<ListboxSelect>>", self.select_set)

        set_new_button = tk.Button(
            self.page, width=4, text="New:", command=self.add_set
        )
        set_new_button.grid(row=20, column=10, columnspan=1)

        set_copy_button = tk.Button(
            self.page, width=4, text="Copy:", command=self.copy_set
        )
        set_copy_button.grid(row=20, column=11, columnspan=1)

        self.set_new_entry = tk.Entry(self.page, width=8)
        self.set_new_entry.grid(row=21, column=10, columnspan=2)

        set_delete_button = tk.Button(
            self.page, width=8, text="Delete:", command=self.delete_set
        )
        set_delete_button.grid(row=22, column=10, columnspan=2)

        #########################################################

        keys_label = tk.Label(self.page, text="Keys")
        keys_label.grid(row=7, column=12, columnspan=2)

        self.keys_listbox = tk.Listbox(
            self.page, width=12, height=9, exportselection=False
        )
        self.keys_listbox.grid(row=8, column=12, columnspan=2, rowspan=5)
        self.keys_listbox.bind("<<ListboxSelect>>", self.select_key)

        keys_new_button = tk.Button(
            self.page, width=4, text="New:", command=self.add_key
        )
        keys_new_button.grid(row=20, column=12, columnspan=1)

        keys_copy_button = tk.Button(
            self.page, width=4, text="Copy:", command=self.copy_key
        )
        keys_copy_button.grid(row=20, column=13, columnspan=1)

        self.keys_new_entry = tk.Entry(self.page, width=8)
        self.keys_new_entry.grid(row=21, column=12, columnspan=2)

        keys_delete_button = tk.Button(
            self.page, width=8, text="Delete:", command=self.delete_key
        )
        keys_delete_button.grid(row=22, column=12, columnspan=2)

        #########################################################

        values_label = tk.Label(self.page, text="Values")
        values_label.grid(row=7, column=14, columnspan=2)

        self.values_listbox = tk.Listbox(
            self.page, width=12, height=9, exportselection=False
        )
        self.values_listbox.grid(row=8, column=14, columnspan=2, rowspan=5)

        values_new_button = tk.Button(
            self.page, width=8, text="New:", command=self.add_values
        )
        values_new_button.grid(row=20, column=14, columnspan=2)

        self.values_new_entry = tk.Entry(self.page, width=8)
        self.values_new_entry.grid(row=21, column=14, columnspan=2)

        values_delete_button = tk.Button(
            self.page, width=8, text="Delete:", command=self.delete_value
        )
        values_delete_button.grid(row=22, column=14, columnspan=2)

    def add_set(self):
        new_record = self.set_new_entry.get()

        if not new_record:
            return

        self.keyword_set[new_record] = {}
        self.set_listbox.insert(0, new_record)
        self.set_listbox.selection_clear(0, "end")
        self.set_listbox.selection_set(0)
        self.select_set(0)

        self.set_new_entry.delete(0, "end")

    def copy_set(self):
        new_record = self.set_new_entry.get()

        if not new_record:
            return

        selected_set = self.set_listbox.curselection()
        if not selected_set:
            messagebox.showwarning(title="Error!", message="Please select the set")
            return

        selected_set = self.set_listbox.get(selected_set)

        self.keyword_set[new_record] = deepcopy(self.keyword_set[selected_set])
        self.set_listbox.insert(0, new_record)
        self.set_listbox.selection_clear(0, "end")
        self.set_listbox.selection_set(0)
        self.select_set(0)

        self.set_new_entry.delete(0, "end")

    def add_key(self):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            messagebox.showwarning(title="Error!", message="Please select the set")
            return

        selected_set = self.set_listbox.get(selected_set)

        new_record = self.keys_new_entry.get()

        if not new_record:
            return

        self.keyword_set[selected_set][new_record] = []
        self.keys_listbox.insert(0, new_record)
        self.keys_listbox.selection_clear(0, "end")
        self.keys_listbox.selection_set(0)
        self.select_key(0)

        self.keys_new_entry.delete(0, "end")

    def copy_key(self):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            messagebox.showwarning(title="Error!", message="Please select the set")
            return

        selected_set = self.set_listbox.get(selected_set)

        selected_key = self.keys_listbox.curselection()
        if not selected_key:
            messagebox.showwarning(title="Error!", message="Please select the key")
            return

        selected_key = self.keys_listbox.get(selected_key)

        new_record = self.keys_new_entry.get()

        if not new_record:
            return

        self.keyword_set[selected_set][new_record] = deepcopy(
            self.keyword_set[selected_set][selected_key]
        )
        self.keys_listbox.insert(0, new_record)
        self.keys_listbox.selection_clear(0, "end")
        self.keys_listbox.selection_set(0)
        self.select_key(0)

        self.keys_new_entry.delete(0, "end")

    def add_values(self):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            messagebox.showwarning(title="Error!", message="Please select the set")
            return

        selected_set = self.set_listbox.get(selected_set)

        selected_key = self.keys_listbox.curselection()
        if not selected_key:
            messagebox.showwarning(title="Error!", message="Please select the key")
            return

        selected_key = self.keys_listbox.get(selected_key)

        new_record = self.values_new_entry.get()

        if not new_record:
            return

        if new_record in self.keyword_set[selected_set][selected_key]:
            return

        self.keyword_set[selected_set][selected_key] = [new_record] + self.keyword_set[
            selected_set
        ][selected_key]
        self.values_listbox.insert(0, new_record)
        self.values_listbox.selection_clear(0, "end")
        self.values_listbox.selection_set(0)

        self.values_new_entry.delete(0, "end")

    def select_set(self, arg):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            print("nihuhu")
            return

        selected_set = self.set_listbox.get(selected_set)

        self.keys_listbox.delete(0, "end")
        self.values_listbox.delete(0, "end")

        for key in self.keyword_set[selected_set]:
            self.keys_listbox.insert("end", key)

    def select_key(self, arg):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            return

        selected_set = self.set_listbox.get(selected_set)

        selected_key = self.keys_listbox.curselection()
        if not selected_key:
            print("nihuhu")
            return

        selected_key = self.keys_listbox.get(selected_key)
        self.values_listbox.delete(0, "end")

        for value in self.keyword_set[selected_set][selected_key]:
            self.values_listbox.insert("end", value)

    def delete_set(self):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            print("nihuhu")
            return

        selected_set_name = self.set_listbox.get(selected_set)
        self.set_listbox.delete(selected_set)

        self.keys_listbox.delete(0, "end")
        self.values_listbox.delete(0, "end")

        del self.keyword_set[selected_set_name]

    def delete_key(self):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            return

        selected_set = self.set_listbox.get(selected_set)

        selected_key = self.keys_listbox.curselection()
        if not selected_key:
            print("nihuhu")
            return

        selected_key_name = self.keys_listbox.get(selected_key)
        self.keys_listbox.delete(selected_key)

        self.values_listbox.delete(0, "end")
        del self.keyword_set[selected_set][selected_key_name]

    def delete_value(self):
        selected_set = self.set_listbox.curselection()
        if not selected_set:
            return

        selected_set = self.set_listbox.get(selected_set)

        selected_key = self.keys_listbox.curselection()
        if not selected_key:
            print("nihuhu")
            return

        selected_key = self.keys_listbox.get(selected_key)

        selected_value = self.values_listbox.curselection()
        if not selected_value:
            return

        selected_value_name = self.values_listbox.get(selected_value)
        self.values_listbox.delete(selected_value)

        self.keyword_set[selected_set][selected_key].remove(selected_value_name)

    def insert_keyword_set(self):
        self.set_listbox.delete(0, "end")
        self.keys_listbox.delete(0, "end")
        self.values_listbox.delete(0, "end")

        for set_name in self.keyword_set:
            self.set_listbox.insert("end", set_name)

    def grid(self):
        self.grid_basic_labels()
        self.grid_host()
        self.grid_guest()
        self.grid_complex()
        self.grid_write()
        self.grid_slurm_section()
        self.grid_set_keywords_section()

        self.grid_gaussian_route_section()
        self.grid_gaussian_additional_input()

    def get_state(self):
        state = {
            "inputBegin": "",
            "additionalInput": "",
            "slurmConfig": "",
            "keywordSet": {},
        }

        state["inputBegin"] = self.route_section.get("1.0", "end")
        state["additionalInput"] = self.additional_section.get("1.0", "end")
        state["slurmConfig"] = self.slurm_text_g16.get("1.0", "end")
        state["keywordSet"] = self.keyword_set

        return state

    def load_state(self, state):
        if "inputBegin" in state:
            self.route_section.delete("1.0", "end")
            self.route_section.insert("end", state["inputBegin"])

        if "additionalInput" in state:
            self.additional_section.delete("1.0", "end")
            self.additional_section.insert("end", state["additionalInput"])

        if "slurmConfig" in state:
            self.slurm_text_g16.delete("1.0", "end")
            self.slurm_text_g16.insert("end", state["slurmConfig"])

        if "keywordSet" in state:
            self.keyword_set = state["keywordSet"]
            self.insert_keyword_set()
