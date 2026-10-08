"""
Created on Tue May 21 16:05:03 2019

@author: michal
"""

import json
from os.path import join

try:
    from pymol import cmd
except ImportError:
    pass

try:
    import paramiko
except ImportError:
    pass

import tkinter as tk
from tkinter import ttk
from tkinter import simpledialog
from tkinter import filedialog
from tkinter import messagebox


class JobStatusGUI:
    def __init__(self, page):
        self.page = page

        self.ntbk = ttk.Notebook(self.page, height=700, width=1320)

        self.job_monitor = ttk.Frame(self.ntbk)
        self.login_data = ttk.Frame(self.ntbk)

        self.ntbk.add(self.job_monitor, text="Job status")
        self.ntbk.add(self.login_data, text="Login data")

        self.ntbk.grid(column=0, row=0, columnspan=20)

        self.tree_headers = ["ID", "Path", "Script", "Status", "Time", "Comment"]
        self.tree_headers2width = {
            "ID": 90,
            "Path": 500,
            "Script": 140,
            "Status": 80,
            "Time": 100,
            "Comment": 200,
        }

        self.client = paramiko.client.SSHClient()
        self.client.load_system_host_keys()
        self.client.set_missing_host_key_policy(paramiko.RejectPolicy)

        self.connected = False
        self.current_selection_tree = None

        self.custom_buttons_no = 18
        self.custom_buttons_per_row = 6
        self.custom_buttons = []
        self.custom_buttons_data = []

        self.actual_status = {}

    def grid_job_monitor(self):
        self.tree_data = ttk.Treeview(
            self.job_monitor, columns=self.tree_headers, show="headings", heigh=15
        )
        for header in self.tree_headers:
            self.tree_data.heading(header, text=header)
            self.tree_data.column(header, width=self.tree_headers2width[header])
        self.tree_data.grid(row=0, column=0, columnspan=20, rowspan=15)
        self.tree_data.bind("<Button-1>", self.set_dir)

        column_no = 21
        get_status_button = tk.Button(
            self.job_monitor, text="Get status", width=15, command=self.get_status
        )
        get_status_button.grid(row=0, column=column_no, columnspan=2)

        cancel_job_button = tk.Button(
            self.job_monitor, text="Cancel job", width=15, command=self.scancel
        )
        cancel_job_button.grid(row=1, column=column_no, columnspan=2)

        forget_button = tk.Button(
            self.job_monitor, text="Forget", width=15, command=self.sremove_py
        )
        forget_button.grid(row=2, column=column_no, columnspan=2)

        self.filter_entry = tk.Entry(self.job_monitor, width=7)
        self.filter_entry.grid(row=3, column=column_no)

        filter_button = tk.Button(
            self.job_monitor, width=5, text="*", command=self.filter_jobs
        )
        filter_button.grid(row=3, column=column_no + 1)

        directory_view_label = tk.Label(self.job_monitor, text="Directory contains:")
        directory_view_label.grid(row=20, column=0, columnspan=2)

        self.directory_view_list = tk.Listbox(self.job_monitor, width=40, height=15)
        self.directory_view_list.grid(row=21, column=0, columnspan=2, rowspan=15)

        refresh_button = tk.Button(
            self.job_monitor,
            text="Refresh",
            width=20,
            command=self.refresh_directory_view,
        )
        refresh_button.grid(row=21, column=2)

        download_button = tk.Button(
            self.job_monitor, text="Download", width=20, command=self.download_file
        )
        download_button.grid(row=22, column=2)

        to_pymol_button = tk.Button(
            self.job_monitor,
            text="to Pymol",
            width=20,
            command=self.download_and_load_to_pymol,
        )
        to_pymol_button.grid(row=23, column=2)

        output_label = tk.Label(self.job_monitor, text="Command output")
        output_label.grid(row=20, column=3, columnspan=4)

        self.output_text = tk.Text(self.job_monitor, width=80, height=16)
        self.output_text.grid(row=21, column=3, columnspan=4, rowspan=15)

        row_actual = 45
        col_actual = 0

        for i in range(self.custom_buttons_no):
            new_button = tk.Button(
                self.job_monitor,
                width=20,
                command=lambda arg=i: self.custom_button_command(arg),
            )
            new_button.grid(row=row_actual, column=col_actual)
            new_button.bind(
                "<Button-3>", lambda e, arg=i: self.custom_button_set(e, arg)
            )
            self.custom_buttons.append(new_button)
            self.custom_buttons_data.append({})
            col_actual += 1
            if col_actual >= self.custom_buttons_per_row:
                col_actual = 0
                row_actual += 1

    def custom_button_command(self, button_ind):
        if not self.connected:
            messagebox.showwarning(
                title="Cannot execute",
                message="You have to connect to host before command execution",
            )
            return

        if "command" not in self.custom_buttons_data[button_ind]:
            messagebox.showwarning(
                title="Cannot execute", message="No command for this button"
            )
            return

        if not self.custom_buttons_data[button_ind]["command"]:
            messagebox.showwarning(
                title="Cannot execute", message="No command for this button"
            )
            return

        command2execute = self.custom_buttons_data[button_ind]["command"]
        current_sel = self.tree_data.focus()
        if current_sel == "":
            messagebox.showwarning(title="Cannot execute", message="Please select job")
            return

        dir2go = self.tree_data.item(current_sel)["values"][1]

        file_selection = self.directory_view_list.curselection()

        if not file_selection:
            messagebox.showwarning(
                title="Cannot execute", message="Please select file to execute"
            )
            return

        file_selection = self.directory_view_list.get(file_selection)

        stdin, stdout, stderr = self.client.exec_command(
            "cd " + dir2go + " ; " + command2execute + " " + file_selection
        )

        output = "".join(list(stdout.readlines()))

        self.output_text.delete("1.0", "end")
        self.output_text.insert("end", output)

    def custom_button_set(self, event, button_ind):
        if "text" in self.custom_buttons_data[button_ind]:
            new_button_name = simpledialog.askstring(
                title="Button name",
                prompt="Select button name",
                initialvalue=self.custom_buttons_data[button_ind]["text"],
            )
        else:
            new_button_name = simpledialog.askstring(
                title="Button name", prompt="Select button name"
            )

        if not new_button_name:
            return

        self.custom_buttons[button_ind].config(text=new_button_name)
        self.custom_buttons_data[button_ind]["text"] = new_button_name

        if "command" in self.custom_buttons_data[button_ind]:
            new_button_command = simpledialog.askstring(
                title="Button command",
                prompt="Select button command",
                initialvalue=self.custom_buttons_data[button_ind]["command"],
            )
        else:
            new_button_command = simpledialog.askstring(
                title="Button command", prompt="Select button command"
            )

        self.custom_buttons_data[button_ind]["command"] = new_button_command

    def set_dir(self, event):
        item = self.tree_data.identify_row(event.y)

        if item:
            info = self.tree_data.item(item, "values")
            self.current_selection_tree = info

            if self.connected:
                dir2print = info[1]

                stdin, stdout, stderr = self.client.exec_command("ls -p " + dir2print)
                files_list = list(stdout.readlines())

                self.directory_view_list.delete(0, "end")
                for filename in files_list:
                    self.directory_view_list.insert("end", filename.strip())

                self.output_text.delete("1.0", "end")

    def get_status(self):
        if not self.connected:
            messagebox.showwarning(
                title="Cannot get status!",
                message="You have to be connected with host to get actual status",
            )

        job_manager_dir = self.job_manager_dir_entry.get()
        if job_manager_dir[-1] != "/":
            job_manager_dir += "/"

        command = " python " + job_manager_dir + "squeuePy.py -json"

        stdin, stdout, stderr = self.client.exec_command(command)

        result = list(stdout.readlines())
        result = " ".join(result)
        result = result.replace("'", '"')
        status = json.loads(result)

        self.actual_status = status

        self.tree_data.delete(*self.tree_data.get_children())
        self.directory_view_list.delete(0, "end")
        self.output_text.delete("1.0", "end")

        for main_key in status:
            result_list = status[main_key]
            for row in result_list:
                table_row = (
                    row["jobID"],
                    row["RunningDir"],
                    row["Script file"],
                    row["Status"],
                    row["Time"],
                    row["Comment"],
                )
                self.tree_data.insert("", "end", values=table_row)

    def filter_jobs(self):
        filter_key = self.filter_entry.get()

        self.tree_data.delete(*self.tree_data.get_children())
        self.directory_view_list.delete(0, "end")
        self.output_text.delete("1.0", "end")

        for main_key in self.actual_status:
            result_list = self.actual_status[main_key]
            for row in result_list:
                string_row = (
                    row["jobID"]
                    + row["RunningDir"]
                    + row["Script file"]
                    + row["Comment"]
                )
                if filter_key in string_row:
                    table_row = (
                        row["jobID"],
                        row["RunningDir"],
                        row["Script file"],
                        row["Status"],
                        row["Time"],
                        row["Comment"],
                    )
                    self.tree_data.insert("", "end", values=table_row)

    def scancel(self):
        if not self.connected:
            messagebox.showwarning(
                title="Cannot scancel!",
                message="You have to be connected with host to cancel job",
            )

        current_sel = self.tree_data.focus()
        if current_sel == "":
            messagebox.showwarning(title="Cannot execute", message="Please select job")
            return

        job_id = self.tree_data.item(current_sel)["values"][0]

        job_manager_dir = self.job_manager_dir_entry.get()
        if job_manager_dir[-1] != "/":
            job_manager_dir += "/"

        command = "scancel " + str(job_id)

        stdin, stdout, stderr = self.client.exec_command(command)

    def sremove_py(self):
        if not self.connected:
            messagebox.showwarning(
                title="Cannot forget!",
                message="You have to be connected with host to forget job",
            )

        current_sel = self.tree_data.focus()
        if current_sel == "":
            messagebox.showwarning(title="Cannot execute", message="Please select job")
            return

        job_id = self.tree_data.item(current_sel)["values"][0]

        job_manager_dir = self.job_manager_dir_entry.get()
        if job_manager_dir[-1] != "/":
            job_manager_dir += "/"

        command = " python " + job_manager_dir + "sremove.py " + str(job_id)

        stdin, stdout, stderr = self.client.exec_command(command)

        item2forget = self.tree_data.selection()[0]
        self.tree_data.delete(item2forget)

    def refresh_directory_view(self):
        if not self.connected:
            messagebox.showwarning(
                title="Cannot execute!", message="You have to be connected with host"
            )

        current_sel = self.tree_data.focus()
        if current_sel == "":
            messagebox.showwarning(title="Cannot execute", message="Please select row")
            return

        item = self.tree_data.focus()
        info = self.tree_data.item(item, "values")
        self.current_selection_tree = info

        if self.connected:
            dir2print = info[1]

            stdin, stdout, stderr = self.client.exec_command("ls -p " + dir2print)
            files_list = list(stdout.readlines())

            self.directory_view_list.delete(0, "end")
            for filename in files_list:
                self.directory_view_list.insert("end", filename.strip())

            self.output_text.delete("1.0", "end")

    def download_file(self):
        if not self.connected:
            messagebox.showwarning(
                title="Cannot execute",
                message="You have to connect to host before command execution",
            )
            return

        current_sel = self.tree_data.focus()
        if current_sel == "":
            messagebox.showwarning(title="Cannot execute", message="Please select job")
            return

        dir2go = self.tree_data.item(current_sel)["values"][1]

        file_selection = self.directory_view_list.curselection()

        if not file_selection:
            messagebox.showwarning(
                title="Cannot execute", message="Please select file to execute"
            )
            return

        file_selection = self.directory_view_list.get(file_selection)

        full_path = join(dir2go, file_selection)

        sftp = self.client.open_sftp()

        sftp.get(full_path, file_selection)

        sftp.close()

    def download_and_load_to_pymol(self):
        if not self.connected:
            messagebox.showwarning(
                title="Cannot execute",
                message="You have to connect to host before command execution",
            )
            return

        current_sel = self.tree_data.focus()
        if current_sel == "":
            messagebox.showwarning(title="Cannot execute", message="Please select job")
            return

        dir2go = self.tree_data.item(current_sel)["values"][1]

        file_selection = self.directory_view_list.curselection()

        if not file_selection:
            messagebox.showwarning(
                title="Cannot execute", message="Please select file to execute"
            )
            return

        file_selection = self.directory_view_list.get(file_selection)

        full_path = join(dir2go, file_selection)

        sftp = self.client.open_sftp()

        sftp.get(full_path, file_selection)

        sftp.close()
        cmd.load(file_selection)

    def grid_login_data(self):
        login_label = tk.Label(self.login_data, text="login")
        login_label.grid(row=0, column=0)

        self.login_entry = tk.Entry(self.login_data, width=20)
        self.login_entry.grid(row=0, column=1)

        host_label = tk.Label(self.login_data, text="host")
        host_label.grid(row=1, column=0)

        self.host_entry = tk.Entry(self.login_data, width=20)
        self.host_entry.grid(row=1, column=1)

        port_label = tk.Label(self.login_data, text="port")
        port_label.grid(row=2, column=0)

        self.port_entry = tk.Entry(self.login_data, width=20)
        self.port_entry.grid(row=2, column=1)
        self.port_entry.insert(0, "22")

        password_label = tk.Label(self.login_data, text="password")
        password_label.grid(row=3, column=0)

        self.password_entry = tk.Entry(self.login_data, width=20)
        self.password_entry.grid(row=3, column=1)

        job_manager_dir_label = tk.Label(self.login_data, text="JobManagerPro dir")
        job_manager_dir_label.grid(row=4, column=0)

        self.job_manager_dir_entry = tk.Entry(self.login_data, width=20)
        self.job_manager_dir_entry.grid(row=4, column=1)

        connect_button = tk.Button(
            self.login_data, width=20, text="Connect", command=self.connect
        )
        connect_button.grid(row=5, column=0)

        disconnect_button = tk.Button(
            self.login_data, width=20, text="Disconnect", command=self.disconnect
        )
        disconnect_button.grid(row=5, column=1)

        status_label = tk.Label(self.login_data, text="Status")
        status_label.grid(row=6, column=0)

        self.status_entry = tk.Entry(self.login_data, width=20)
        self.status_entry.grid(row=6, column=1)

        self.status_entry.insert(0, "Disconnected")
        self.status_entry.configure(state="readonly")

        download_label = tk.Label(self.login_data, text="Download dir")
        download_label.grid(row=7, column=0)

        download_dir_button = tk.Button(
            self.login_data, text="Change", width=20, command=self.change_download_dir
        )
        download_dir_button.grid(row=7, column=1)

        self.download_entry = tk.Entry(self.login_data, width=60)
        self.download_entry.grid(row=8, column=0, columnspan=3)
        self.download_entry.configure(state="readonly")

    def change_download_dir(self):
        new_dir = filedialog.askdirectory()
        if not new_dir:
            return

        self.download_entry.configure(state="normal")
        self.download_entry.delete(0, "end")
        self.download_entry.insert("end", new_dir)
        self.download_entry.configure(state="readonly")

    def connect(self):
        host = self.host_entry.get()
        login = self.login_entry.get()
        port = int(self.port_entry.get())
        password = self.password_entry.get()

        try:
            self.client.connect(host, port=port, username=login, password=password)

            self.status_entry.configure(state="normal")
            self.status_entry.delete(0, "end")
            self.status_entry.insert(0, "Connected")
            self.status_entry.configure(state="readonly")
            self.connected = True
        except Exception:
            messagebox.showwarning(
                title="Connection error!",
                message="Cannot connect to host! Please check login, password "
                "and internet connection",
            )

    def disconnect(self):
        if self.connected:
            self.client.close()

            self.status_entry.configure(state="normal")
            self.status_entry.delete(0, "end")
            self.status_entry.insert(0, "Disconnected")
            self.status_entry.configure(state="readonly")
            self.connected = False

    def get_state(self):
        state = {}

        state["login"] = self.login_entry.get()
        state["host"] = self.host_entry.get()
        state["port"] = self.port_entry.get()
        state["jobManagerDir"] = self.job_manager_dir_entry.get()
        state["password"] = self.password_entry.get()
        state["customButtons"] = self.custom_buttons_data
        state["downloadDir"] = self.download_entry.get()

        return state

    def load_state(self, state):
        self.login_entry.delete(0, "end")
        self.login_entry.insert(0, state["login"])

        self.host_entry.delete(0, "end")
        self.host_entry.insert(0, state["host"])

        self.port_entry.delete(0, "end")
        self.port_entry.insert(0, state["port"])

        self.password_entry.delete(0, "end")
        self.password_entry.insert(0, state["password"])

        self.job_manager_dir_entry.delete(0, "end")
        self.job_manager_dir_entry.insert(0, state["jobManagerDir"])

        if "customButtons" in state:
            self.custom_buttons_data = state["customButtons"]
            self.refresh_custom_buttons()

        if "downloadDir" in state:
            self.download_entry.configure(state="normal")
            self.download_entry.delete(0, "end")
            self.download_entry.insert("end", state["downloadDir"])
            self.download_entry.configure(state="readonly")

    def refresh_custom_buttons(self):
        for i, data in enumerate(self.custom_buttons_data):
            if "text" in data:
                self.custom_buttons[i].config(text=data["text"])

        len_diff = len(self.custom_buttons) - len(self.custom_buttons_data)
        if len_diff > 0:
            for i in range(len_diff):
                self.custom_buttons_data.append({})

    def grid(self):
        self.grid_job_monitor()
        self.grid_login_data()
