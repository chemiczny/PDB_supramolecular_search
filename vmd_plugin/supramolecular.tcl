# GUI for testing the supramolecular search
package require Tk
package provide supramolecular 0.0

namespace eval ::Supramolecular:: {
	variable window 
	variable lastSessionFile [ file join [ file dirname [ file normalize [ info script ] ] ] .lastSession ]
	variable xyzDir ""
	variable logFile ""
	variable ligand_anions_anion2file
	variable ligand_ligand2file
	variable aminoacid_anion2file
	variable aminoacid_aminoacid2file
	variable currentXYZ ""
	variable lastLoadedStructureId -1
	variable currentMoleculeContext
	variable lastMoleculeContext

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
	  # $window.tabs add [ frame $window.tabs.pdbAnalyse ] -text "PDB analyser"
	  $window.tabs add [ frame $window.tabs.logAnalyse ] -text "XYZ analyser"

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

	  listbox $window.tabs.logAnalyse.list_structures -height 20 -width 100
	  grid $window.tabs.logAnalyse.list_structures -row 3 -column 3 -rowspan 3
	  bind $window.tabs.logAnalyse.list_structures <<ListboxSelect>> { ::Supramolecular::printStructure }
	  # scrollbar $window.tabs.logAnalyse.sb -command [list $window.tabs.logAnalyse.list_structures yview] 
	  # $window.tabs.logAnalyse.list_structures configure -yscrollcommand [list $window.tabs.logAnalyse.sb set]
	  # grid $window.tabs.logAnalyse.sb -row 3 -column 4 -rowspan 3

	  button $window.tabs.logAnalyse.bShowPiacid -text "Show only Pi-acid" -command showPiAcid
	  grid $window.tabs.logAnalyse.bShowPiacid -row 0 -column 1

	  button $window.tabs.logAnalyse.bShowPiacidAndAnion -text "Show Pi-acid and anion" -command showPiAcidAndAnion
	  grid $window.tabs.logAnalyse.bShowPiacidAndAnion -row 0 -column 2

	  button $window.tabs.logAnalyse.bShowFull -text "Show full structure" -command showFullStructure
	  grid $window.tabs.logAnalyse.bShowFull -row 1 -column 1

	  button $window.tabs.logAnalyse.bShowPiacidAndEnv -text "Show Pi-acid and env" -command showPiAcidAndEnv
	  grid $window.tabs.logAnalyse.bShowPiacidAndEnv -row 1 -column 2

	  button $window.tabs.logAnalyse.bShowInteractions -text "Show interactions" -command ::Supramolecular::showInteractions
	  grid $window.tabs.logAnalyse.bShowInteractions -row 0 -column 3

	  if { [ file exists $lastSessionFile ] == 1  } {
	   	  readLastSession $lastSessionFile

	   	  if { $xyzDir != "" } {
	   	  	puts $xyzDir
	   	  	fillLists $xyzDir
	   	  }
	   }
	}

	proc loadLogFile { } {
	  variable logFile
	  set logTypes { 
	  	{ {Log files} {.log} }
	  	{ {Out files} {.out} }
	  	{ {Txt files} {.txt} }
	  	{ {All files} * }  
	  }
	  
	  set filename [ tk_getOpenFile -filetypes $logTypes ]

	  if { $filename ne "" } {
	  	set logFile $filename
	  	updateLastSessionFile logFile $filename
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
		updateLastSessionFile xyzDir $dir
	}

	proc updateLastSessionFile  { parameterName parameterValue } {
		variable lastSessionFile

		set sessionF [ open $lastSessionFile r ]
		set sessionData [ list ]

		while { [gets $sessionF data] >= 0 } {
			lappend sessionData $data
		}

		close $sessionF

		set sessionF [ open $lastSessionFile w ]
		set newParameterFound -1

		foreach line $sessionData {
			if { [ string match "* $parameterName *" $line ]} {
				set newParameterFound 1
				puts $sessionF "set $parameterName $parameterValue"
			} else {
				puts $sessionF $line
			}
		}
		if { $newParameterFound < 0  } {
			puts $sessionF "set $parameterName $parameterValue"
		}

		close $sessionF

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
			set ligand_ligand2file [ getDictFromDir [ file join $dir "ligands_ENV"  ]  ]
			insertSortedListIntoWidget $ligand_ligand2file $window.tabs.logAnalyse.list_ligands
			set aminoacid_anion2file [ getDictFromDir [ file join $dir "aminoacids_anions"  ]  ]
			insertSortedListIntoWidget $aminoacid_anion2file $window.tabs.logAnalyse.list_aminoacids_anions
			set aminoacid_aminoacid2file [ getDictFromDir [ file join $dir "aminoacids_ENV"  ]  ]
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
		$widget delete 0 end
		foreach ligand [ lsort [ dict keys $inputDict ] ] {
			$widget insert end $ligand
		}
	}

	proc readLastSession { lastSessionFile } {
		variable xyzDir
		variable logFile
		source $lastSessionFile
	}

	proc handleStructures { structureList } {
		variable window
		variable ligand_anions_anion2file
		variable ligand_ligand2file
		variable aminoacid_anion2file
		variable aminoacid_aminoacid2file
		variable currentXYZ
		variable currentMoleculeContext

		set xyzKey [ $window.tabs.logAnalyse.$structureList curselection ]
		puts $xyzKey
		set currentDict ""

		if { [ llength $xyzKey ] > 1 || $xyzKey == ""   } {
			return
		}

		set currentMoleculeContext $structureList
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
			if { [ string match "*Model*" $data ] } {
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
			if { [ string match "*Model*" $data ] } {
				if { $commentsNo == $structureNo } {
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

		mol delrep 0 top 
		mol representation Licorice
		mol selection "not name X"
		mol addrep top

		mol material Steel
		mol representation Lines
		mol selection "name X"
		mol addrep top
	}

	proc getCurrentStructureInfo { } {
		set currF [ open temp.xyz r ]
		gets $currF data 
		gets $currF data
		close $currF

		return $data
	}

	proc showFullStructure { } {

	}

	proc showPiAcid { } {

	}

	proc showPiAcidAndEnv { } {
		variable currentMoleculeContext



	}

	proc showPiAcidAndAnion { } {

	}

	proc showInteractions { } {
		set structureInfo [ getCurrentStructureInfo ]
		set structureInfo [ xyzComment2dict $structureInfo ]

		findInteractionInLog $structureInfo
	}

	proc xyzComment2dict { comment } {
		return [ string map { : "" _ " " }  $comment ]
	}

	proc findInteractionInLog { structureInfo } {
		variable logFile
		puts $logFile
		puts $structureInfo
		
		set pdbCode [ dict get $structureInfo PDBCode ]
		set modelNo [ dict get $structureInfo ModelNo ]
		set piAcidId [ dict get $structureInfo Pi-acid-id ]
		set anionId -1
		if { [ dict exists $structureInfo Anion-id ] } {
			set anionId [ dict get $structureInfo Anion-id ]
		}

		set log [ open $logFile r ]

		if { $anionId > 0 } {
			while { [ gets $log data ] >= 0 } {
				set pdbCoderead [ lindex $data 0 ]
				set modelNoread [ lindex $data 19 ]
				set piAcidIdread [ lindex $data 3 ]
				set anionIdread [ lindex $data 6 ]
				if { $pdbCoderead == $pdbCode && $modelNoread == $modelNo && $piAcidIdread == $piAcidId && $anionIdread == $anionId } {
					printInteraction $data
				}
			}
		} else {
			while { [ gets $log data ] >= 0 } {
				set pdbCoderead [ lindex $data 0 ]
				set modelNoread [ lindex $data 19 ]
				set piAcidIdread [ lindex $data 3 ]

				if { $pdbCoderead == $pdbCode && $modelNoread == $modelNo && $piAcidIdread == $piAcidId } {
					printInteraction $data
				}
			}
		}
		
		close $log
	}
}

proc printInteraction { data } {
	set centroid [ lrange $data 13 15 ]
	set anionAtom [ lrange $data 16 18 ]
	set anionName [ lindex $data 4 ]

	set anionTextPoint [ list ]

	foreach c $centroid a $anionAtom {
		lappend anionTextPoint [ expr $a + 0.18*( $c - $a )  ]
	}

	draw line $centroid $anionAtom 
	draw text $anionTextPoint $anionName

	set middlePoint [ list ]
	foreach c $centroid a $anionAtom {
		lappend middlePoint [ expr ( $c + $a )/2. ]
	}

	set distance  [ string range [lindex $data 9 ] 0 3 ] 
	draw text $middlePoint $distance
}


::Supramolecular::supramolecular_gui