"""
Created on Sun Jul 22 19:40:59 2018

@author: michal
"""

try:
    from pymol import plugins
except ImportError:
    pass

import tkinter as tk
from tkinter import filedialog
from tkinter import ttk

from supramolecular_composition import SupramolecularComposition


def fetch_dialog(simulation=False):
    if simulation:
        root = tk.Tk()
    else:
        root = plugins.get_tk_root()

    app_data = {"cifDir": False}

    main_window = tk.Toplevel(root)
    main_window.title("Supramolecular analyser")
    main_window.minsize(1335, 400)

    main_window.resizable(1, 1)

    canvas = tk.Canvas(main_window, width=1320, height=900)
    canvas.grid(row=0, column=0, columnspan=50)

    main_scrollbar = tk.Scrollbar(main_window, orient="vertical", command=canvas.yview)
    main_scrollbar.grid(row=0, column=50, rowspan=1)

    def move_down(event):
        canvas.yview_scroll(1, "units")

    def move_up(event):
        canvas.yview_scroll(-1, "units")

    canvas.configure(yscrollcommand=main_scrollbar.set)
    canvas.configure(scrollregion=(0, 0, 1320, 1800))
    canvas.bind_all("<Down>", move_down)
    canvas.bind_all("<Up>", move_up)

    main_frame = tk.Frame(canvas, width=1320, height=900)
    canvas.create_window((0, 0), window=main_frame, anchor="nw")

    nb = ttk.Notebook(main_frame, height=700, width=1320)

    page_anion_pi = ttk.Frame(nb)
    page_pi_pi = ttk.Frame(nb)
    page_cation_pi = ttk.Frame(nb)
    page_anion_cation = ttk.Frame(nb)
    page_h_bonds = ttk.Frame(nb)
    page_metal_ligand = ttk.Frame(nb)
    page_linear_anion_pi = ttk.Frame(nb)
    page_planar_anion_pi = ttk.Frame(nb)
    page_qm = ttk.Frame(nb)
    page_job_status = ttk.Frame(nb)

    nb.add(page_anion_pi, text="AnionPi")
    nb.add(page_pi_pi, text="PiPi")
    nb.add(page_cation_pi, text="CationPi")
    nb.add(page_anion_cation, text="AnionCation")
    nb.add(page_h_bonds, text="HBonds")
    nb.add(page_metal_ligand, text="MetalLigand")
    nb.add(page_linear_anion_pi, text="LinearAnionPi")
    nb.add(page_planar_anion_pi, text="PlanarAnionPi")
    nb.add(page_qm, text="QM assistant")
    nb.add(page_job_status, text="Job monitor")

    nb.grid(column=0, row=0, columnspan=20)

    supramolecular_composition = SupramolecularComposition(
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
    )

    ######################
    # GENERAL
    ######################

    def read_log_dir():
        app_data["logDir"] = filedialog.askdirectory()
        if not app_data["logDir"]:
            print("nihuhu")
            return

        ent_log_dir.configure(state="normal")
        ent_log_dir.delete(0, "end")
        ent_log_dir.insert(0, app_data["logDir"])
        ent_log_dir.configure(state="readonly")

        supramolecular_composition.read_all_logs_from_dir(app_data["logDir"])

    but_read_log_dir = tk.Button(
        main_frame, width=20, text="Read log Dir", command=read_log_dir
    )
    but_read_log_dir.grid(row=2, column=0, columnspan=2)

    ent_log_dir = tk.Entry(main_frame, width=45)
    ent_log_dir.configure(state="readonly")
    ent_log_dir.grid(row=2, column=2, columnspan=4)

    def select_cif():
        app_data["cifDir"] = filedialog.askdirectory()
        if not app_data["cifDir"]:
            print("nihuhu")
            return

        ent_cif_dir.configure(state="normal")
        ent_cif_dir.delete(0, "end")
        ent_cif_dir.insert(0, app_data["cifDir"])
        ent_cif_dir.configure(state="readonly")

        supramolecular_composition.select_cif_dir(app_data["cifDir"])

    but_cif_dir = tk.Button(main_frame, width=10, command=select_cif, text="Cif dir")
    but_cif_dir.grid(row=2, column=6)

    ent_cif_dir = tk.Entry(main_frame, width=45)
    ent_cif_dir.configure(state="readonly")
    ent_cif_dir.grid(row=2, column=7, columnspan=3)

    def set_parallel_selection():
        state = var_parallel_selection.get()
        supramolecular_composition.set_parallel_selection(state)

    lab_parallel_selection = tk.Label(main_frame, text="Parallel selection")
    lab_parallel_selection.grid(row=2, column=10)

    var_parallel_selection = tk.IntVar()
    chk_parallel_selection = tk.Checkbutton(
        main_frame, variable=var_parallel_selection, command=set_parallel_selection
    )
    chk_parallel_selection.grid(row=2, column=11)

    action_menu = {}

    column = 1

    lab_use_page = tk.Label(main_frame, width=10, text="Use:")
    lab_use_page.grid(row=4, column=0)

    lab_do_not_use_page = tk.Label(main_frame, width=10, text="Exclude:")
    lab_do_not_use_page.grid(row=5, column=0)

    for label in supramolecular_composition.action_labels:
        action_menu[label] = {}

        action_menu[label]["label"] = tk.Label(main_frame, text=label)
        action_menu[label]["label"].grid(row=3, column=column)

        actual_row = 4

        action_menu[label]["checkValue"] = tk.IntVar()
        action_menu[label]["checkbox"] = tk.Checkbutton(
            main_frame, variable=action_menu[label]["checkValue"]
        )
        action_menu[label]["checkbox"].grid(row=actual_row, column=column)

        actual_row += 1

        action_menu[label]["checkValueExclude"] = tk.IntVar()
        action_menu[label]["checkboxExclude"] = tk.Checkbutton(
            main_frame, variable=action_menu[label]["checkValueExclude"]
        )
        action_menu[label]["checkboxExclude"].grid(row=actual_row, column=column)

        column += 1

    lab_show_int = tk.Label(main_frame, width=10, text="Show:")
    lab_show_int.grid(row=6, column=0)

    show_menu = {}
    column = 1

    for label in supramolecular_composition.action_labels:
        show_menu[label] = {}

        show_menu[label]["checkValue"] = tk.IntVar()
        show_menu[label]["checkbox"] = tk.Checkbutton(
            main_frame, variable=show_menu[label]["checkValue"]
        )
        show_menu[label]["checkbox"].grid(row=6, column=column)

        column += 1

    def merge_results():
        supramolecular_composition.merge(action_menu)

    ent_recording_state = tk.Entry(main_frame, width=20)
    ent_recording_state.grid(row=3, column=9, columnspan=2)
    ent_recording_state.insert(0, "Not recording")
    ent_recording_state.configure(state="readonly")

    def start_recording():
        ent_recording_state.configure(state="normal")
        ent_recording_state.delete(0, "end")
        ent_recording_state.insert(0, "Recording")
        ent_recording_state.configure(state="readonly")

        supramolecular_composition.start_recording()

    but_start_recording = tk.Button(
        main_frame, width=20, text="Start recording", command=start_recording
    )
    but_start_recording.grid(row=3, column=11, columnspan=2)

    def stop_recording():
        ent_recording_state.configure(state="normal")
        ent_recording_state.delete(0, "end")
        ent_recording_state.insert(0, "Not recording")
        ent_recording_state.configure(state="readonly")

        supramolecular_composition.stop_recording()

    but_stop_recording = tk.Button(
        main_frame, width=20, text="Stop recording", command=stop_recording
    )
    but_stop_recording.grid(row=4, column=11, columnspan=2)

    but_merge = tk.Button(main_frame, width=20, text="Merge!", command=merge_results)
    but_merge.grid(row=5, column=9, columnspan=2)

    def show_all_interactions():
        supramolecular_composition.show_all(show_menu)

    but_show_many = tk.Button(
        main_frame, width=20, text="Show", command=show_all_interactions
    )
    but_show_many.grid(row=6, column=9, columnspan=2)

    but_save_state = tk.Button(
        main_frame,
        width=20,
        text="Save GUI state",
        command=supramolecular_composition.save_state,
    )
    but_save_state.grid(row=5, column=11, columnspan=2)

    but_load_state = tk.Button(
        main_frame,
        width=20,
        text="Load GUI state",
        command=supramolecular_composition.load_state,
    )
    but_load_state.grid(row=6, column=11, columnspan=2)

    # INDIVIDUALS
    lab_anion_groups_as_individuals = tk.Label(
        main_frame, width=30, text="Anion's groups as indyviduals:"
    )
    lab_anion_groups_as_individuals.grid(row=7, column=1, columnspan=3)

    chk_anion_groups_as_individuals = tk.Checkbutton(
        main_frame, variable=supramolecular_composition.anions_as_individuals
    )
    chk_anion_groups_as_individuals.grid(row=7, column=4)

    lab_rings_as_individuals = tk.Label(
        main_frame, width=30, text="Residues's rings as indyviduals:"
    )
    lab_rings_as_individuals.grid(row=7, column=5, columnspan=2)

    chk_rings_as_individuals = tk.Checkbutton(
        main_frame, variable=supramolecular_composition.residue_rings_as_individuals
    )
    chk_rings_as_individuals.grid(row=7, column=7)
    ######################
    # ALL
    ######################

    supramolecular_composition.grid()

    if simulation:
        main_frame.mainloop()
