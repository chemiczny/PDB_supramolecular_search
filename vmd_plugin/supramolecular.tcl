# GUI for testing the supramolecular search
package require Tk
package provide supramolecular 0.0

namespace eval ::Supramolecular:: {
	variable window 
	variable lastSessionFile [ file join [ file dirname [ file normalize [ info script ] ] ] .lastSession ]
	variable xyzDir ""
	variable ligand_anions_anion2file
	variable ligand_ligand2file
	variable aminoacid_anion2file
	variable aminoacid_aminoacid2file
	variable currentXYZ ""
	variable lastLoadedStructureId -1

	proc supramolecular_gui { } {
		variable window
		variable lastSessionFile
		variable xyzDir

		if { [winfo exist .supramolecular] } {
	      wm deiconify $window
	      raise $window
	      return
	   }

	   set window [toplevel .supremolecular]
	   wm title $window "Supramolecular helper"

	  ttk::notebook $window.tabs
	  pack $window.tabs
	  $window.tabs add [ frame $window.tabs.pdbAnalyse ] -text "PDB analyser"
	  $window.tabs add [ frame $window.tabs.logAnalyse ] -text "Log analyser"

	  button $window.tabs.logAnalyse.bLog -text "Load log file" -command ::Supramolecular::loadLogFile
	  button $window.tabs.logAnalyse.bXyz -text "Load xyz directory" -command ::Supramolecular::loadXyzDir

	  grid $window.tabs.logAnalyse.bLog -row 0 -column 0
	  grid $window.tabs.logAnalyse.bXyz -row 1 -column 0

	  label $window.tabs.logAnalyse.lab_ligand -text "Ligands workspace:"
	  grid  $window.tabs.logAnalyse.lab_ligand -row 2 -column 1 -columnspan 2

	  listbox $window.tabs.logAnalyse.list_ligand_anions -height 10
	  grid $window.tabs.logAnalyse.list_ligand_anions -row 3 -column 0
	  bind $window.tabs.logAnalyse.list_ligand_anions <<ListboxSelect>> { ::Supramolecular::handleStructures list_ligand_anions }

	  listbox $window.tabs.logAnalyse.list_ligands -height 10
	  grid $window.tabs.logAnalyse.list_ligands -row 3 -column 1
	  bind $window.tabs.logAnalyse.list_ligands <<ListboxSelect>> { ::Supramolecular::handleStructures list_ligands }

	  label $window.tabs.logAnalyse.lab_aminoacids -text "Aminoacids workspace:"
	  grid  $window.tabs.logAnalyse.lab_aminoacids -row 4 -column 1 -columnspan 2

	  listbox $window.tabs.logAnalyse.list_aminoacids_anions -height 10
	  grid $window.tabs.logAnalyse.list_aminoacids_anions -row 5 -column 0
	  bind $window.tabs.logAnalyse.list_aminoacids_anions <<ListboxSelect>> { ::Supramolecular::handleStructures list_aminoacids_anions }

	  listbox $window.tabs.logAnalyse.list_aminoacids -height 10
	  grid $window.tabs.logAnalyse.list_aminoacids -row 5 -column 1
	  bind $window.tabs.logAnalyse.list_aminoacids <<ListboxSelect>> { ::Supramolecular::handleStructures list_aminoacids }

	  text $window.tabs.logAnalyse.text_currXyz -height 1 -width 40
	  grid $window.tabs.logAnalyse.text_currXyz -row 2 -column 3

	  listbox $window.tabs.logAnalyse.list_structures -height 20 -width 60
	  grid $window.tabs.logAnalyse.list_structures -row 3 -column 3 -rowspan 3
	  bind $window.tabs.logAnalyse.list_structures <<ListboxSelect>> { ::Supramolecular::printStructure }
	  # scrollbar $window.tabs.logAnalyse.sb -command [list $window.tabs.logAnalyse.list_structures yview] 
	  # $window.tabs.logAnalyse.list_structures configure -yscrollcommand [list $window.tabs.logAnalyse.sb set]
	  # grid $window.tabs.logAnalyse.sb -row 3 -column 4 -rowspan 3

	  if { [ file exists $lastSessionFile ] == 1  } {
	   	  readLastSession $lastSessionFile

	   	  if { $xyzDir != "" } {
	   	  	puts $xyzDir
	   	  	fillLists $xyzDir
	   	  }
	   }
	}

	proc loadLogFile { } {
	  set logTypes { 
	  	{ {Log files} {.log} }
	  	{ {Out files} {.out} }
	  	{ {Txt files} {.txt} }
	  	{ {All files} * }  
	  }
	  
	  set filename [ tk_getOpenFile -filetypes $logTypes ]

	  if { $filename ne "" } {
	  	puts lol
	  }
	}

	proc loadXyzDir { } {
		variable xyzDir

		if { $xyzDir != "" } {
			set dir [ tk_chooseDirectory -initialdir $xyzDir ]
		} else {
			set dir [ tk_chooseDirectory ]
		}
	
		fillLists $dir

	}

	proc updateLastSessionFile  { parameterName parameterValue } {
		variable lastSessionFile


	}

	proc fillLists { dir } {
		variable window
		variable ligand_anions_anion2file
		variable ligand_ligand2file
		variable aminoacid_anion2file
		variable aminoacid_aminoacid2file

		if { $dir ne "" } {
			set ligand_anions_anion2file [ getDictFromDir [ file join $dir "ligands_anions"  ]  ]
			insertSortedListIntoWidget $ligand_anions_anion2file $window.tabs.logAnalyse.list_ligand_anions
			set ligand_ligand2file [ getDictFromDir [ file join $dir "ligands"  ]  ]
			insertSortedListIntoWidget $ligand_ligand2file $window.tabs.logAnalyse.list_ligands
			set aminoacid_anion2file [ getDictFromDir [ file join $dir "aminoacids_anions"  ]  ]
			insertSortedListIntoWidget $aminoacid_anion2file $window.tabs.logAnalyse.list_aminoacids_anions
			set aminoacid_aminoacid2file [ getDictFromDir [ file join $dir "aminoacids"  ]  ]
			insertSortedListIntoWidget $aminoacid_aminoacid2file $window.tabs.logAnalyse.list_aminoacids
			
		}
	}

	proc getDictFromDir { dir } {
		puts $dir
		set ligand_anions [ glob -nocomplain -directory $dir *.xyz ]
		set ligand_anions_anion2file [ dict create ]
		foreach ligand_anions_file $ligand_anions {
			set key [ lindex [ split [ file tail $ligand_anions_file ] "_" ] 0 ]
			dict set ligand_anions_anion2file $key $ligand_anions_file
		}
		return $ligand_anions_anion2file
	}

	proc insertSortedListIntoWidget { inputDict widget  } {
		variable window
		foreach ligand [ lsort [ dict keys $inputDict ] ] {
			$widget insert end $ligand
		}
	}

	proc readLastSession { lastSessionFile } {
		variable xyzDir
		source $lastSessionFile
	}

	proc handleStructures { structureList } {
		variable window
		variable ligand_anions_anion2file
		variable ligand_ligand2file
		variable aminoacid_anion2file
		variable aminoacid_aminoacid2file
		variable currentXYZ

		set xyzKey [ $window.tabs.logAnalyse.$structureList curselection ]
		puts $xyzKey
		set currentDict ""

		if { [ llength $xyzKey ] > 1 || $xyzKey == ""   } {
			return
		}

		switch $structureList {
			list_ligand_anions {
				set currentDict $ligand_anions_anion2file
			}
			list_ligands {
				set currentDict $ligand_ligand2file
			}
			list_aminoacids_anions {
				set currentDict $aminoacid_anion2file
			}
			list_aminoacids {
				set currentDict $aminoacid_aminoacid2file
			}
			default {}
		}

		set xyzKey [ lindex [ lsort [ dict keys $currentDict ] ] $xyzKey ]
		set xyzFile [ dict get $currentDict $xyzKey ]
		set currentXYZ $xyzFile
		$window.tabs.logAnalyse.text_currXyz delete 0.0 end
		$window.tabs.logAnalyse.text_currXyz insert 0.0  "$structureList  [ file tail $xyzFile ] "

		set structures [ getDataFromXYZ $xyzFile ]
		$window.tabs.logAnalyse.list_structures delete 0 end
		foreach  structure  $structures {
			$window.tabs.logAnalyse.list_structures insert end $structure
		}
	}

	proc getDataFromXYZ { xyzFile } {
		set fxyz [ open $xyzFile r ]
		set comments [ list ]

		while { [gets $fxyz data] >= 0 } {
			if { [ string match "*Structure*" $data ] } {
				lappend comments $data
			}
		}

		close $fxyz
		return $comments
	}

	proc printStructure { } {
		variable currentXYZ
		variable window
		variable lastLoadedStructureId

		set structureNo [ $window.tabs.logAnalyse.list_structures curselection ]
		if { $structureNo == "" } {
			return
		}

		set fxyz [ open $currentXYZ r ]
		set commentsNo 0
		set previousLine ""
		set actualComment ""

		while { [gets $fxyz data] >= 0 } {
			if { [ string match "*Structure*" $data ] } {
				if { $commentsNo == $structureNo } {
					puts "mam ja!"
					set actualComment $data
					break
				}
				incr commentsNo
			}
			set previousLine $data
		}

		set tempXYZ [ open "temp.xyz" w ]
		puts $tempXYZ $previousLine
		puts $tempXYZ $actualComment
		for {set i 0} {$i < $previousLine} {incr i} {
			gets $fxyz data
			puts $tempXYZ $data
		}

		close $tempXYZ
		close $fxyz

		if { $lastLoadedStructureId >= 0 } {
			mol delete $lastLoadedStructureId
		}
		set lastLoadedStructureId [ mol new temp.xyz ]
		mol representation Licorice
		mol addrep top
	}

}


::Supramolecular::supramolecular_gui