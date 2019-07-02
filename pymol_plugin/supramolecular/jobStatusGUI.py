#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 16:05:03 2019

@author: michal
"""

import sys
import json
from os.path import join
try:
    from pymol import cmd
except:
    pass

try:
    import paramiko
except:
    pass

if sys.version_info[0] < 3:
    import Tkinter
    import tkMessageBox, tkFileDialog
    import tkSimpleDialog as simpledialog
    import ttk
    
else:
    import tkinter as Tkinter
    import tkinter.ttk as ttk
    from tkinter import simpledialog
    from tkinter import filedialog as tkFileDialog
    from tkinter import messagebox as tkMessageBox
    
class JobStatusGUI:
    def __init__(self, page):
        self.page = page
        
        self.ntbk = ttk.Notebook(self.page, height = 700, width = 1320)
        
        self.jobMonitor = ttk.Frame(self.ntbk)
        self.loginData = ttk.Frame(self.ntbk)
        
        self.ntbk.add( self.jobMonitor, text = "Job status" )
        self.ntbk.add( self.loginData, text = "Login data")
        
        self.ntbk.grid(column = 0, row = 0, columnspan = 20)
        
        self.treeHeaders = ["ID", "Path", "Script", "Status", "Time", "Comment"]
        self.treeHeaders2width = {"ID" : 90, "Path" : 500, "Script" : 140, "Status" : 80, "Time" : 100, "Comment" : 200 }
        
        self.client = paramiko.client.SSHClient()
        self.client.load_system_host_keys()
        self.client.set_missing_host_key_policy(paramiko.RejectPolicy)
        
        self.connected = False
        self.currentSelectionTree = None
        
        self.customButtonsNo = 18
        self.customButtonsPerRow = 6
        self.customButtons = []
        self.customButtonsData = []
        
        self.actualStatus = {}
        
    def gridJobMonitor(self):
        self.tree_data = ttk.Treeview(self.jobMonitor, columns = self.treeHeaders , show = "headings", heigh = 15 )
        for header in self.treeHeaders:
            self.tree_data.heading(header, text = header)
            self.tree_data.column(header, width = self.treeHeaders2width[header])
        self.tree_data.grid(row = 0, column = 0, columnspan = 20, rowspan = 15)
        self.tree_data.bind("<Button-1>", self.setDir)
        
        columnNo = 21
        getStatusButton = Tkinter.Button(self.jobMonitor, text = "Get status", width = 15, command = self.getStatus)
        getStatusButton.grid(row=0, column = columnNo, columnspan = 2)
        
        cancelJobButton = Tkinter.Button(self.jobMonitor, text = "Cancel job", width = 15, command = self.scancel)
        cancelJobButton.grid(row=1, column = columnNo, columnspan = 2)
        
        forgetButton = Tkinter.Button(self.jobMonitor, text = "Forget", width = 15, command = self.sremovePy)
        forgetButton.grid(row=2, column = columnNo, columnspan = 2)
        
        self.filterEntry = Tkinter.Entry(self.jobMonitor, width = 7)
        self.filterEntry.grid(row= 3, column = columnNo)
        
        filterButton = Tkinter.Button(self.jobMonitor, width = 5, text = "*", command = self.filterJobs)
        filterButton.grid(row = 3, column = columnNo + 1)
        
        directoryViewLabel = Tkinter.Label(self.jobMonitor, text = "Directory contains:")
        directoryViewLabel.grid( row = 20, column = 0 , columnspan = 2)
        
        self.directoryViewList = Tkinter.Listbox(self.jobMonitor, width = 40, height = 15 )
        self.directoryViewList.grid(row = 21, column = 0, columnspan = 2, rowspan = 15)
        
        refreshButton = Tkinter.Button(self.jobMonitor, text = "Refresh", width = 20, command = self.refreshDirectoryView )
        refreshButton.grid(row= 21, column =2)
        
        downloadButton = Tkinter.Button(self.jobMonitor, text = "Download", width = 20, command = self.downloadFile)
        downloadButton.grid(row=22, column = 2)
        
        toPymolButton = Tkinter.Button(self.jobMonitor, text = "to Pymol", width = 20, command = self.downloadAndLoadToPymol )
        toPymolButton.grid(row = 23, column = 2)
        
        outputLabel = Tkinter.Label(self.jobMonitor, text = "Command output")
        outputLabel.grid(row = 20, column = 3, columnspan = 4)
        
        self.outputText = Tkinter.Text(self.jobMonitor, width = 80, height =  16)
        self.outputText.grid(row = 21, column = 3, columnspan = 4, rowspan = 15)
        
        rowActual = 45
        colActual = 0
        
        for i in range( self.customButtonsNo):
            newButton = Tkinter.Button( self.jobMonitor, width = 20, command = lambda arg = i : self.customButtonCommand(arg) )
            newButton.grid(row = rowActual, column = colActual)
            newButton.bind("<Button-3>", lambda e, arg = i:self.customButtonSet(e, arg))
            self.customButtons.append(newButton)
            self.customButtonsData.append({})
            colActual += 1
            if colActual >= self.customButtonsPerRow:
                colActual = 0
                rowActual += 1
        
    def customButtonCommand(self, buttonInd):
        if not self.connected:
            tkMessageBox.showwarning(title = "Cannot execute", message = "You have to connect to host before command execution")
            return
        
        if not "command" in self.customButtonsData[buttonInd]:
            tkMessageBox.showwarning(title = "Cannot execute", message = "No command for this button")
            return
        
        if not self.customButtonsData[buttonInd]["command"]:
            tkMessageBox.showwarning(title = "Cannot execute", message = "No command for this button")
            return
        
        command2execute = self.customButtonsData[buttonInd]["command"]
        currentSel = self.tree_data.focus()
        if currentSel == "" :
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select job")
            return
        
        dir2go = self.tree_data.item(currentSel)["values"][1]
        
        fileSelection = self.directoryViewList.curselection()
        
        if not fileSelection:
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select file to execute")
            return
        
        fileSelection = self.directoryViewList.get(fileSelection)
        
        stdin, stdout, stderr = self.client.exec_command("cd "+dir2go + " ; " + command2execute + " "+ fileSelection)
        
        output = "".join(list(stdout.readlines()))
        
        self.outputText.delete("1.0", "end")
        self.outputText.insert("end", output)
        
    def customButtonSet(self, event, buttonInd):
        if "text" in self.customButtonsData[buttonInd]:
            newButtonName = simpledialog.askstring(title = "Button name", prompt = "Select button name", initialvalue = self.customButtonsData[buttonInd]["text"])
        else:
            newButtonName = simpledialog.askstring(title = "Button name", prompt = "Select button name")
        
        if not newButtonName:
            return
        
        self.customButtons[buttonInd].config(text = newButtonName)
        self.customButtonsData[buttonInd]["text"] = newButtonName
        
        if "command" in self.customButtonsData[buttonInd]:
            newButtonCommand = simpledialog.askstring(title = "Button command", prompt = "Select button command", initialvalue = self.customButtonsData[buttonInd]["command"])
        else:
            newButtonCommand = simpledialog.askstring(title = "Button command", prompt = "Select button command")
        
        self.customButtonsData[buttonInd]["command"] = newButtonCommand
        
    def setDir(self, event):
        item = self.tree_data.identify_row(event.y)
        
        if item:
            info = self.tree_data.item(item, 'values')
            self.currentSelectionTree = info
            
            if self.connected:
                dir2print = info[1]
                
                stdin, stdout, stderr = self.client.exec_command("ls -p "+dir2print)
                filesList = list(stdout.readlines())
                
                self.directoryViewList.delete(0, "end")
                for filename in filesList:
                    self.directoryViewList.insert("end", filename.strip())
                    
                self.outputText.delete("1.0", "end")
                
        
    def getStatus(self):
        if not self.connected:
            tkMessageBox.showwarning(title = "Cannot get status!", message = "You have to be connected with host to get actual status")
            
        jobManagerDir = self.jobManagerDirEntry.get()
        if jobManagerDir[-1] != "/":
            jobManagerDir += "/"
        
        command = " python " + jobManagerDir + "squeuePy.py -json"
        
        stdin, stdout, stderr = self.client.exec_command(command)
        
        result = list(stdout.readlines())
        result = " ".join(result)
        result = result.replace("'", '"')
        status = json.loads( result )
        
        self.actualStatus = status
        
        self.tree_data.delete(*self.tree_data.get_children())
        self.directoryViewList.delete(0, "end")
        self.outputText.delete("1.0", "end")
        
        for mainKey in status:
            resultList = status[mainKey]
            for row in resultList:
                tableRow = ( row["jobID"] , row["RunningDir"], row["Script file"], row["Status"], row["Time"], row["Comment"] )
                self.tree_data.insert('', "end" , values = tableRow  )
                
    def filterJobs(self):
        filterKey = self.filterEntry.get()
        
        self.tree_data.delete(*self.tree_data.get_children())
        self.directoryViewList.delete(0, "end")
        self.outputText.delete("1.0", "end")
        
        for mainKey in self.actualStatus:
            resultList = self.actualStatus[mainKey]
            for row in resultList:
                stringRow = row["jobID"] + row["RunningDir"] + row["Script file"] + row["Comment"]
                if filterKey in stringRow:
                    tableRow = ( row["jobID"] , row["RunningDir"], row["Script file"], row["Status"], row["Time"], row["Comment"] )
                    self.tree_data.insert('', "end" , values = tableRow  )
                
    def scancel(self):
        if not self.connected:
            tkMessageBox.showwarning(title = "Cannot scancel!", message = "You have to be connected with host to cancel job")
            
        currentSel = self.tree_data.focus()
        if currentSel == "" :
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select job")
            return
        
        jobID = self.tree_data.item(currentSel)["values"][0]
            
        jobManagerDir = self.jobManagerDirEntry.get()
        if jobManagerDir[-1] != "/":
            jobManagerDir += "/"
        
        command = "scancel " + str(jobID)
        
        stdin, stdout, stderr = self.client.exec_command(command)
    
    def sremovePy(self):
        if not self.connected:
            tkMessageBox.showwarning(title = "Cannot forget!", message = "You have to be connected with host to forget job")
            
        currentSel = self.tree_data.focus()
        if currentSel == "" :
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select job")
            return
        
        jobID = self.tree_data.item(currentSel)["values"][0]
            
        jobManagerDir = self.jobManagerDirEntry.get()
        if jobManagerDir[-1] != "/":
            jobManagerDir += "/"
        
        command =  " python " + jobManagerDir + "sremove.py " + str(jobID)
        
        stdin, stdout, stderr = self.client.exec_command(command)
        
        item2forget = self.tree_data.selection()[0]
        self.tree_data.delete(item2forget)
    
    def refreshDirectoryView(self):
        if not self.connected:
            tkMessageBox.showwarning(title = "Cannot execute!", message = "You have to be connected with host")
            
        currentSel = self.tree_data.focus()
        if currentSel == "" :
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select row")
            return
        
        item = self.tree_data.focus()
        info = self.tree_data.item(item, 'values')
        self.currentSelectionTree = info
        
        if self.connected:
            dir2print = info[1]
            
            stdin, stdout, stderr = self.client.exec_command("ls -p "+dir2print)
            filesList = list(stdout.readlines())
            
            self.directoryViewList.delete(0, "end")
            for filename in filesList:
                self.directoryViewList.insert("end", filename.strip())
                
            self.outputText.delete("1.0", "end")
    
    def downloadFile(self):
        if not self.connected:
            tkMessageBox.showwarning(title = "Cannot execute", message = "You have to connect to host before command execution")
            return
        
        currentSel = self.tree_data.focus()
        if currentSel == "" :
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select job")
            return
        
        dir2go = self.tree_data.item(currentSel)["values"][1]
        
        fileSelection = self.directoryViewList.curselection()
        
        if not fileSelection:
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select file to execute")
            return
        
        fileSelection = self.directoryViewList.get(fileSelection)
        
        fullPath = join(dir2go, fileSelection)
        
        sftp = self.client.open_sftp()
        
        sftp.get(fullPath, fileSelection)
        
        sftp.close()
    
    def downloadAndLoadToPymol(self):
        if not self.connected:
            tkMessageBox.showwarning(title = "Cannot execute", message = "You have to connect to host before command execution")
            return
        
        currentSel = self.tree_data.focus()
        if currentSel == "" :
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select job")
            return
        
        dir2go = self.tree_data.item(currentSel)["values"][1]
        
        fileSelection = self.directoryViewList.curselection()
        
        if not fileSelection:
            tkMessageBox.showwarning(title = "Cannot execute", message = "Please select file to execute")
            return
        
        fileSelection = self.directoryViewList.get(fileSelection)
        
        fullPath = join(dir2go, fileSelection)
        
        sftp = self.client.open_sftp()
        
        sftp.get(fullPath, fileSelection)
        
        sftp.close()
        cmd.load(fileSelection)
    
    def gridLoginData(self):
        loginLabel = Tkinter.Label(self.loginData, text = "login")
        loginLabel.grid(row = 0, column = 0)
        
        self.loginEntry = Tkinter.Entry(self.loginData, width = 20)
        self.loginEntry.grid(row = 0, column = 1)
        
        hostLabel = Tkinter.Label(self.loginData, text = "host")
        hostLabel.grid(row= 1, column = 0)
        
        self.hostEntry = Tkinter.Entry(self.loginData, width = 20)
        self.hostEntry.grid(row = 1, column = 1)
        
        portLabel = Tkinter.Label(self.loginData, text = "port")
        portLabel.grid(row = 2, column = 0)
        
        self.portEntry = Tkinter.Entry(self.loginData, width = 20)
        self.portEntry.grid(row = 2, column = 1)
        self.portEntry.insert(0, "22")
        
        passwordLabel = Tkinter.Label(self.loginData, text = "password")
        passwordLabel.grid(row = 3, column =0)
        
        self.passwordEntry = Tkinter.Entry(self.loginData, width = 20)
        self.passwordEntry.grid(row = 3, column = 1)
        
        jobManagerDirLabel = Tkinter.Label(self.loginData, text = "JobManagerPro dir")
        jobManagerDirLabel.grid(row = 4, column  = 0)
        
        self.jobManagerDirEntry = Tkinter.Entry(self.loginData, width = 20)
        self.jobManagerDirEntry.grid( row = 4, column = 1 )
        
        connectButton = Tkinter.Button(self.loginData, width = 20, text = "Connect", command = self.connect)
        connectButton.grid(row = 5, column = 0)
        
        disconnectButton = Tkinter.Button(self.loginData, width = 20, text = "Disconnect", command = self.disconnect)
        disconnectButton.grid(row = 5, column = 1)
        
        statusLabel = Tkinter.Label(self.loginData, text = "Status")
        statusLabel.grid(row = 6, column = 0)
        
        self.statusEntry = Tkinter.Entry(self.loginData, width = 20)
        self.statusEntry.grid(row = 6, column = 1)
        
        self.statusEntry.insert(0, "Disconnected")
        self.statusEntry.configure(state = "readonly")
        
        downloadLabel = Tkinter.Label(self.loginData, text = "Download dir")
        downloadLabel.grid(row = 7, column = 0)
        
        downloadDirButton = Tkinter.Button(self.loginData, text = "Change", width = 20, command = self.changeDownloadDir)
        downloadDirButton.grid(row = 7, column = 1)
        
        self.downloadEntry = Tkinter.Entry(self.loginData, width = 60)
        self.downloadEntry.grid(row = 8, column = 0, columnspan = 3)
        self.downloadEntry.configure(state = "readonly")
        
    def changeDownloadDir(self):
        newDir =  tkFileDialog.askdirectory()
        if not newDir:
            return
        
        self.downloadEntry.configure(state = "normal")
        self.downloadEntry.delete(0, "end")
        self.downloadEntry.insert("end", newDir)
        self.downloadEntry.configure(state = "readonly")
    
    def connect(self):
        host = self.hostEntry.get()
        login = self.loginEntry.get()
        port =  int(self.portEntry.get())
        password = self.passwordEntry.get()
        
        try:
            self.client.connect(host, port=port , username=login, password=password)
            
            self.statusEntry.configure(state = "normal")
            self.statusEntry.delete(0, "end")
            self.statusEntry.insert(0, "Connected")
            self.statusEntry.configure(state = "readonly")
            self.connected = True
        except:
            tkMessageBox.showwarning(title = "Connection error!", message = "Cannot connect to host! Please check login, password and internet connection")
    
    def disconnect(self):
        if self.connected:
            self.client.close()
            
            self.statusEntry.configure(state = "normal")
            self.statusEntry.delete(0, "end")
            self.statusEntry.insert(0, "Disconnected")
            self.statusEntry.configure(state = "readonly")
            self.connected = False
    
    def getState(self):
        state = {}
        
        state["login"] = self.loginEntry.get()
        state["host"] = self.hostEntry.get()
        state["port"] = self.portEntry.get()
        state["jobManagerDir"] = self.jobManagerDirEntry.get()
        state["password"] = self.passwordEntry.get()
        state["customButtons"] = self.customButtonsData
        state["downloadDir"] = self.downloadEntry.get()
        
        return state
    
    def loadState(self, state):
        self.loginEntry.delete(0, "end")
        self.loginEntry.insert(0, state["login"])
        
        self.hostEntry.delete(0, "end")
        self.hostEntry.insert(0, state["host"])
        
        self.portEntry.delete(0, "end")
        self.portEntry.insert(0, state["port"])
        
        self.passwordEntry.delete(0, "end")
        self.passwordEntry.insert(0, state["password"])
        
        self.jobManagerDirEntry.delete(0, "end")
        self.jobManagerDirEntry.insert(0, state["jobManagerDir"])
        
        if "customButtons" in state:
            self.customButtonsData = state["customButtons"]
            self.refreshCustomButtons()
            
        if "downloadDir" in state:
            self.downloadEntry.configure("normal")
            self.downloadEntry.delete(0, "end")
            self.downloadEntry.insert("end", state["downloadDir"])
            self.downloadEntry.configure(state = "readonly")

    def refreshCustomButtons(self):
        for i, data in enumerate(self.customButtonsData):
            if "text" in data:
                self.customButtons[i].config(text = data["text"])

    def grid(self):
        self.gridJobMonitor()
        self.gridLoginData()