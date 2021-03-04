proc ConstraintVectorTable {FileVar TableId TableDict CondName VarName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $TableId MyTableId
        upvar $TableDict MyTableDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set AuxList [list]
            if {[lindex [lindex $Groups $i] 6] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME ${VarName}_X"
                set Table [lindex [lindex $Groups $i] 7]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
            }
            if {[lindex [lindex $Groups $i] 11] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME ${VarName}_Y"
                set Table [lindex [lindex $Groups $i] 12]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 0
            }
            if {[lindex [lindex $Groups $i] 16] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table2 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME ${VarName}_Z"
                set Table [lindex [lindex $Groups $i] 17]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table2 0
            }
            dict set MyTableDict [lindex [lindex $Groups $i] 1] TableList $AuxList
        }
    }
}

# TODO: it may be dangerous to write Tables without format (puts -nonewline $FileVar [format  "%.10f" [lindex $Table $j]])

#-------------------------------------------------------------------------------

proc PressureTable {FileVar TableId TableDict CondName VarName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $TableId MyTableId
        upvar $TableDict MyTableDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set AuxList [list]
            if {[lindex [lindex $Groups $i] 3] eq "Phreatic_Line"} {
                if {[lindex [lindex $Groups $i] 22] eq "Table_Interpolation"} {
                    incr MyTableId
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                    lappend AuxList $MyTableId
                    puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                    set Table [lindex [lindex $Groups $i] 23]
                    for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                        puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                    }
                    puts $MyFileVar "End Table"
                    puts $MyFileVar ""
                } else {
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
                }
                if {[lindex [lindex $Groups $i] 24] eq "Table_Interpolation"} {
                    incr MyTableId
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 $MyTableId
                    lappend AuxList $MyTableId
                    puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                    set Table [lindex [lindex $Groups $i] 25]
                    for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                        puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                    }
                    puts $MyFileVar "End Table"
                    puts $MyFileVar ""
                } else {
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 0
                }

            } elseif {[lindex [lindex $Groups $i] 3] eq "Phreatic_Surface"} {
                if {[lindex [lindex $Groups $i] 22] eq "Table_Interpolation"} {
                    incr MyTableId
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                    lappend AuxList $MyTableId
                    puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                    set Table [lindex [lindex $Groups $i] 23]
                    for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                        puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                    }
                    puts $MyFileVar "End Table"
                    puts $MyFileVar ""
                } else {
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
                }
                if {[lindex [lindex $Groups $i] 24] eq "Table_Interpolation"} {
                    incr MyTableId
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 $MyTableId
                    lappend AuxList $MyTableId
                    puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                    set Table [lindex [lindex $Groups $i] 25]
                    for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                        puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                    }
                    puts $MyFileVar "End Table"
                    puts $MyFileVar ""
                } else {
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 0
                }
                if {[lindex [lindex $Groups $i] 26] eq "Table_Interpolation"} {
                    incr MyTableId
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table2 $MyTableId
                    lappend AuxList $MyTableId
                    puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                    set Table [lindex [lindex $Groups $i] 27]
                    for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                        puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                    }
                    puts $MyFileVar "End Table"
                    puts $MyFileVar ""
                } else {
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table2 0
                }
            } else {
                if {[lindex [lindex $Groups $i] 22] eq "Table_Interpolation"} {
                    incr MyTableId
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                    lappend AuxList $MyTableId
                    puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                    set Table [lindex [lindex $Groups $i] 23]
                    for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                        puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                    }
                    puts $MyFileVar "End Table"
                    puts $MyFileVar ""
                } else {
                    dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
                }
            }
            dict set MyTableDict [lindex [lindex $Groups $i] 1] TableList $AuxList
        }
    }
}

#-------------------------------------------------------------------------------

proc VectorTable {FileVar TableId TableDict CondName VarName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $TableId MyTableId
        upvar $TableDict MyTableDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set AuxList [list]
            if {[lindex [lindex $Groups $i] 5] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME ${VarName}_X"
                set Table [lindex [lindex $Groups $i] 6]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
            }
            if {[lindex [lindex $Groups $i] 9] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME ${VarName}_Y"
                set Table [lindex [lindex $Groups $i] 10]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 0
            }
            if {[lindex [lindex $Groups $i] 13] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table2 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME ${VarName}_Z"
                set Table [lindex [lindex $Groups $i] 14]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table2 0
            }
            dict set MyTableDict [lindex [lindex $Groups $i] 1] TableList $AuxList
        }
    }
}

#-------------------------------------------------------------------------------

proc NormalTangentialTable {FileVar TableId TableDict CondName NormalVarName TangentialVarName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $TableId MyTableId
        upvar $TableDict MyTableDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set AuxList [list]
            if {[lindex [lindex $Groups $i] 9] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME $NormalVarName"
                set Table [lindex [lindex $Groups $i] 10]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
            }
            if {[lindex [lindex $Groups $i] 13] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME $TangentialVarName"
                set Table [lindex [lindex $Groups $i] 14]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table1 0
            }
            dict set MyTableDict [lindex [lindex $Groups $i] 1] TableList $AuxList
        }
    }
}

#-------------------------------------------------------------------------------

proc ScalarTable {FileVar TableId TableDict CondName VarName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $TableId MyTableId
        upvar $TableDict MyTableDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set AuxList [list]
            if {[lindex [lindex $Groups $i] 4] eq "Table_Interpolation"} {
                incr MyTableId
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 $MyTableId
                lappend AuxList $MyTableId
                puts $MyFileVar "Begin Table $MyTableId TIME $VarName"
                set Table [lindex [lindex $Groups $i] 5]
                for {set j 2} {$j <= [lindex $Table 1]} {incr j 2} {
                    puts $MyFileVar "  [lindex $Table $j] [lindex $Table [expr { $j+1 }]]"
                }
                puts $MyFileVar "End Table"
                puts $MyFileVar ""
            } else {
                dict set MyTableDict [lindex [lindex $Groups $i] 1] Table0 0
            }
            dict set MyTableDict [lindex [lindex $Groups $i] 1] TableList $AuxList
        }
    }
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------

proc WriteAnchorElements {FileVar Group ElemType ElemName PropertyId ConnectivityType} {
    set Entities [GiD_EntitiesGroups get [lindex $Group 1] elements -element_type $ElemType]
    if {[llength $Entities] > 0} {
        upvar $FileVar MyFileVar
        #puts $MyFileVar $Entities
        set Nodes [GiD_Info Mesh Nodes]
        set MaxDist -1e20
        set MinDist 1e20
        
        set HighValue 1e10
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
            #puts $MyFileVar [lindex $Entities $j]
            set ElementInfo [GiD_Mesh get element [lindex $Entities $j]]
            set NodeNrOne [lindex $ElementInfo 3]
            set NodeNrTwo [lindex $ElementInfo 4]
            
            set XCoordN1 [lindex $Nodes [expr { ($NodeNrOne-1)*4 + 1}]]
            set YCoordN1 [lindex $Nodes [expr { ($NodeNrOne-1)*4 + 2}]]
            set ZCoordN1 [lindex $Nodes [expr { ($NodeNrOne-1)*4 + 3}]]
                    
            set XCoordN2 [lindex $Nodes [expr { ($NodeNrTwo-1)*4 + 1}]]
            set YCoordN2 [lindex $Nodes [expr { ($NodeNrTwo-1)*4 + 2}]]
            set ZCoordN2 [lindex $Nodes [expr { ($NodeNrTwo-1)*4 + 3}]]
            
            set CurrentDistN1 [expr { sqrt(($XCoordN1 - $HighValue)*($XCoordN1 - $HighValue) + ($YCoordN1 - $HighValue)*($YCoordN1 - $HighValue) + ($ZCoordN1 - $HighValue)*($ZCoordN1 - $HighValue)) }]
            set CurrentDistN2 [expr { sqrt(($XCoordN2 - $HighValue)*($XCoordN2 - $HighValue) + ($YCoordN2 - $HighValue)*($YCoordN2 - $HighValue) + ($ZCoordN2 - $HighValue)*($ZCoordN2 - $HighValue)) }]
            
            if {$CurrentDistN1 > $MaxDist} { 
                set MaxDist $CurrentDistN1
                set MaxNodeNr $NodeNrOne
            }
            if {$CurrentDistN1 < $MinDist} {
                set MinDist $CurrentDistN1
                set MinNodeNr $NodeNrOne
            }
            
            if {$CurrentDistN2 > $MaxDist} { 
                set MaxDist $CurrentDistN2
                set MaxNodeNr $NodeNrTwo
            }
            
            if {$CurrentDistN2 < $MinDist} {
                set MinDist $CurrentDistN2
                set MinNodeNr $NodeNrTwo
            }
        }
        puts $MyFileVar "Begin Elements $ElemName"
        puts $MyFileVar "  [lindex $Entities 0]  $PropertyId  $MinNodeNr $MaxNodeNr"
        # for {set j 0} {$j < [llength $Entities]} {incr j} {
            # puts $MyFileVar "  [lindex $Entities $j]  $PropertyId  [$ConnectivityType [lindex $Entities $j]]"
        # }
        puts $MyFileVar "End Elements"
        puts $MyFileVar ""
        
        
        #set Groups [GiD_Info conditions $CondName groups]
        #puts [lindex $Group 1]  nodes
        #GiD_EntitiesGroups assign [lindex $Group 1] [lindex $Entities 0]
        #GiD_EntitiesGroups assign [lindex $Group 1] lines [lindex $Entities 0]]
        
                    # puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # # Tables
            # puts $MyFileVar "  Begin SubModelPartTables"
            # puts $MyFileVar "  End SubModelPartTables"
            # # Nodes
            # set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            # puts $MyFileVar "  Begin SubModelPartNodes"
            # for {set j 0} {$j < [llength $Entities]} {incr j} {
                # puts $MyFileVar "    [lindex $Entities $j]"
            #}
        
    }
    return [list [lindex $Entities 0] $PropertyId  $MinNodeNr $MaxNodeNr $Entities]

}


proc WriteElements {FileVar Group ElemType ElemName PropertyId ConnectivityType} {
    set Entities [GiD_EntitiesGroups get [lindex $Group 1] elements -element_type $ElemType]
    if {[llength $Entities] > 0} {
        upvar $FileVar MyFileVar

        puts $MyFileVar "Begin Elements $ElemName"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            puts $MyFileVar "  [lindex $Entities $j]  $PropertyId  [$ConnectivityType [lindex $Entities $j]]"
        }
        puts $MyFileVar "End Elements"
        puts $MyFileVar ""
    }
}

#-------------------------------------------------------------------------------

proc WritePropUnionElements {FileVar PropertyId} {
    upvar $FileVar MyFileVar

    set Groups [GiD_Groups list PropagationUnion_3d_6]
    set ElementId [GiD_Info Mesh MaxNumElements]
    set PropUnionElementList [list]

    puts $MyFileVar "Begin Elements UPwSmallStrainInterfaceElement3D6N"
    for {set i 0} {$i < [llength $Groups]} {incr i 6} {
        incr ElementId
        lappend PropUnionElementList $ElementId
        set Connectivities [list [lindex [GiD_EntitiesGroups get "PropagationUnion_3d_6//[lindex $Groups $i]" nodes] 0] \
                                 [lindex [GiD_EntitiesGroups get "PropagationUnion_3d_6//[lindex $Groups [expr {($i+1)}]]" nodes] 0] \
                                 [lindex [GiD_EntitiesGroups get "PropagationUnion_3d_6//[lindex $Groups [expr {($i+2)}]]" nodes] 0] \
                                 [lindex [GiD_EntitiesGroups get "PropagationUnion_3d_6//[lindex $Groups [expr {($i+3)}]]" nodes] 0] \
                                 [lindex [GiD_EntitiesGroups get "PropagationUnion_3d_6//[lindex $Groups [expr {($i+4)}]]" nodes] 0] \
                                 [lindex [GiD_EntitiesGroups get "PropagationUnion_3d_6//[lindex $Groups [expr {($i+5)}]]" nodes] 0]]

        puts $MyFileVar "  $ElementId  $PropertyId  $Connectivities"
    }
    puts $MyFileVar "End Elements"
    puts $MyFileVar ""

    return $PropUnionElementList
}

#-------------------------------------------------------------------------------

proc WriteNodalConditions {FileVar ConditionId ConditionDict Groups CondName PropertyId} {
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionDict MyConditionDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "Begin Conditions $CondName"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                incr MyConditionId
                lappend MyConditionList $MyConditionId
               
                puts $MyFileVar "  $MyConditionId  $PropertyId  [lindex $Entities $j]"
            }
            puts $MyFileVar "End Conditions"
            puts $MyFileVar ""
            dict set MyConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteLineConditions {FileVar ConditionId ConditionDict Groups CondName PropertyDict} {
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionDict MyConditionDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
            puts $MyFileVar "Begin Conditions $CondName"
            for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                incr MyConditionId
                lappend MyConditionList $MyConditionId
                set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                for {set k 0} {$k < [llength $ElementGroup]} {incr k} {
                    if {[dict exists $PropertyDict [lindex $ElementGroup $k]] eq 1} {
                        set PropertyId [dict get $PropertyDict [lindex $ElementGroup $k]]
                        break
                    }
                }
                set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                puts $MyFileVar "  $MyConditionId  $PropertyId  $Connectivities"
            }
            puts $MyFileVar "End Conditions"
            puts $MyFileVar ""
            dict set MyConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    }
}



#-------------------------------------------------------------------------------

proc WriteFaceConditions {FileVar ConditionId ConditionDict Groups CondName PropertyDict} {
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionDict MyConditionDict

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            set MyConditionList [list]
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] faces]
            puts $MyFileVar "Begin Conditions $CondName"
            for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
                incr MyConditionId
                lappend MyConditionList $MyConditionId
                set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
                for {set k 0} {$k < [llength $ElementGroup]} {incr k} {
                    if {[dict exists $PropertyDict [lindex $ElementGroup $k]] eq 1} {
                        set PropertyId [dict get $PropertyDict [lindex $ElementGroup $k]]
                        break
                    }
                }
                set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
                puts $MyFileVar "  $MyConditionId  $PropertyId  $Connectivities"
            }
            puts $MyFileVar "End Conditions"
            puts $MyFileVar ""
            dict set MyConditionDict [lindex [lindex $Groups $i] 1] $MyConditionList
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteTypeFaceConditions {FileVar ConditionId ConditionList Group ElemType CondName PropertyDict} {
    set Entities [GiD_EntitiesGroups get [lindex $Group 1] faces -element_type $ElemType]
    if {[llength [lindex $Entities 1]] > 0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionList MyConditionList

        puts $MyFileVar "Begin Conditions $CondName"
        for {set j 0} {$j < [llength [lindex $Entities 1]]} {incr j} {
            incr MyConditionId
            lappend MyConditionList $MyConditionId
            set ElementGroup [GiD_EntitiesGroups entity_groups element [lindex [lindex $Entities 0] $j]]
            for {set k 0} {$k < [llength $ElementGroup]} {incr k} {
                if {[dict exists $PropertyDict [lindex $ElementGroup $k]] eq 1} {
                    set PropertyId [dict get $PropertyDict [lindex $ElementGroup $k]]
                    break
                }
            }
            set Connectivities [GiD_Mesh get element [lindex [lindex $Entities 0] $j] face [lindex [lindex $Entities 1] $j]]
            puts $MyFileVar "  $MyConditionId  $PropertyId  $Connectivities"
        }
        puts $MyFileVar "End Conditions"
        puts $MyFileVar ""
    }
}


#-------------------------------------------------------------------------------

proc WriteInterfaceConditions {FileVar ConditionId ConditionList Group ElemType CondName PropertyId ConnectivityType} {
    set Entities [GiD_EntitiesGroups get [lindex $Group 1] elements -element_type $ElemType]
    if {[llength $Entities] > 0} {
        upvar $FileVar MyFileVar
        upvar $ConditionId MyConditionId
        upvar $ConditionList MyConditionList

        puts $MyFileVar "Begin Conditions $CondName"
        for {set j 0} {$j < [llength $Entities]} {incr j} {
            incr MyConditionId
            lappend MyConditionList $MyConditionId
            puts $MyFileVar "  $MyConditionId  $PropertyId  [$ConnectivityType [lindex $Entities $j]]"
        }
        puts $MyFileVar "End Conditions"
        puts $MyFileVar ""
    }
}

#-------------------------------------------------------------------------------
proc Line2D2Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4]"
}

#-------------------------------------------------------------------------------
proc Line2D3Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
}

#-------------------------------------------------------------------------------

proc Triangle2D3Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
}

#-------------------------------------------------------------------------------

proc Quadrilateral3D4Connectivities { ElemId } {

    #Note: It is the same for the Tethrahedron3D4

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6]"
}

#-------------------------------------------------------------------------------

proc Quadrilateral2D4Connectivities { ElemId } {

    #Note: It is the same for the Tethrahedron3D4

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6]"
}

#-------------------------------------------------------------------------------

proc Triangle3D3Connectivities { ElemId } {

    #It is the same for the Prism3D6

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
}


#-------------------------------------------------------------------------------

proc Triangle2D6Connectivities { ElemId } {

    #It is the same for the Prism3D6

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]"
}

#-------------------------------------------------------------------------------

proc Hexahedron3D8Connectivities { ElemId } {

    #It is the same for Quadrilateral2D8

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10]"
}

#-------------------------------------------------------------------------------

proc Quadrilateral2D9Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]"
}

#-------------------------------------------------------------------------------

proc Tetrahedron3D10Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]\
    [lindex $ElementInfo 12]"
}

#-------------------------------------------------------------------------------

proc Hexahedron3D20Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]\
    [lindex $ElementInfo 12] [lindex $ElementInfo 13] [lindex $ElementInfo 14]\
    [lindex $ElementInfo 15] [lindex $ElementInfo 16] [lindex $ElementInfo 17]\
    [lindex $ElementInfo 18] [lindex $ElementInfo 19] [lindex $ElementInfo 20]\
    [lindex $ElementInfo 21] [lindex $ElementInfo 22]"
}

#-------------------------------------------------------------------------------

proc Hexahedron3D27Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10] [lindex $ElementInfo 11]\
    [lindex $ElementInfo 12] [lindex $ElementInfo 13] [lindex $ElementInfo 14]\
    [lindex $ElementInfo 15] [lindex $ElementInfo 16] [lindex $ElementInfo 17]\
    [lindex $ElementInfo 18] [lindex $ElementInfo 19] [lindex $ElementInfo 20]\
    [lindex $ElementInfo 21] [lindex $ElementInfo 22] [lindex $ElementInfo 23]\
    [lindex $ElementInfo 24] [lindex $ElementInfo 25] [lindex $ElementInfo 26]\
    [lindex $ElementInfo 27] [lindex $ElementInfo 28] [lindex $ElementInfo 29]"
}

#-------------------------------------------------------------------------------

proc TriangleInterface2D4Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 4]\
    [lindex $ElementInfo 5]"

    #~ set N1(Id) [lindex $ElementInfo 3]
    #~ set N2(Id) [lindex $ElementInfo 4]
    #~ set N3(Id) $N2(Id)
    #~ set N4(Id) [lindex $ElementInfo 5]

    #~ # Obtaining nodes coordinates
    #~ set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    #~ set N1(x) [lindex $NCoord 0]
    #~ set N1(y) [lindex $NCoord 1]
    #~ set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    #~ set N2(x) [lindex $NCoord 0]
    #~ set N2(y) [lindex $NCoord 1]
    #~ #set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    #~ set N3(x) $N2(x)
    #~ set N3(y) $N2(y)
    #~ set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    #~ set N4(x) [lindex $NCoord 0]
    #~ set N4(y) [lindex $NCoord 1]

    #~ # Computing element lengths
    #~ set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 ) }]
    #~ set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 ) }]

    #~ if {$ly <= $lx} {
        #~ set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    #~ } else {
        #~ set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    #~ }

    #~ return $Connectivities
}

#-------------------------------------------------------------------------------

proc QuadrilateralInterface2D4Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6]"

    #~ set N1(Id) [lindex $ElementInfo 3]
    #~ set N2(Id) [lindex $ElementInfo 4]
    #~ set N3(Id) [lindex $ElementInfo 5]
    #~ set N4(Id) [lindex $ElementInfo 6]

    #~ # Obtaining nodes coordinates
    #~ set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    #~ set N1(x) [lindex $NCoord 0]
    #~ set N1(y) [lindex $NCoord 1]
    #~ set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    #~ set N2(x) [lindex $NCoord 0]
    #~ set N2(y) [lindex $NCoord 1]
    #~ set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    #~ set N3(x) [lindex $NCoord 0]
    #~ set N3(y) [lindex $NCoord 1]
    #~ set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    #~ set N4(x) [lindex $NCoord 0]
    #~ set N4(y) [lindex $NCoord 1]

    #~ # Computing element lengths
    #~ set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 ) }]
    #~ set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 ) }]

    #~ if {$ly <= $lx} {
        #~ set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    #~ } else {
        #~ set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    #~ }

    #~ return $Connectivities
}

#-------------------------------------------------------------------------------

proc TetrahedronInterface3D6Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 4] [lindex $ElementInfo 5]"
}

#-------------------------------------------------------------------------------

proc PrismInterface3D6Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]"

    #~ ## Check element orientation (very slow!)
    #~ # Obtaining element volume
    #~ set Volume [lindex [GiD_Info list_entities -more Elements $ElemId] 20]
    #~ set Volume [string trimleft $Volume "Volume="]

    #~ # Check Connectivities
    #~ if {$Volume > -1.0e-10} {
        #~ return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
        #~ [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]"
    #~ } else {
        #~ #W "Reordering nodes of interface element"
        #~ return "[lindex $ElementInfo 3] [lindex $ElementInfo 5] [lindex $ElementInfo 4]\
        #~ [lindex $ElementInfo 6] [lindex $ElementInfo 8] [lindex $ElementInfo 7]"
    #~ }
}

#-------------------------------------------------------------------------------

proc HexahedronInterface3D8Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
    [lindex $ElementInfo 9] [lindex $ElementInfo 10]"

    #~ ## Check element orientation (very slow!)
    #~ # Obtaining element volume
    #~ set Volume [lindex [GiD_Info list_entities -more Elements $ElemId] 22]
    #~ set Volume [string trimleft $Volume "Volume="]

    #~ # Check Connectivities
    #~ if {$Volume > -1.0e-10} {
        #~ return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
        #~ [lindex $ElementInfo 6] [lindex $ElementInfo 7] [lindex $ElementInfo 8]\
        #~ [lindex $ElementInfo 9] [lindex $ElementInfo 10]"
    #~ } else {
        #~ #W "Reordering nodes of interface element"
        #~ return "[lindex $ElementInfo 3] [lindex $ElementInfo 6] [lindex $ElementInfo 5]\
        #~ [lindex $ElementInfo 4] [lindex $ElementInfo 7] [lindex $ElementInfo 10]\
        #~ [lindex $ElementInfo 9] [lindex $ElementInfo 8]"
    #~ }
}

#-------------------------------------------------------------------------------

#proc HexahedronInterface3D8Connectivities { ElemId } {

    ## Obtaining element nodes
    #set ElementInfo [GiD_Mesh get element $ElemId]
    ##ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    #set N1(Id) [lindex $ElementInfo 3]
    #set N2(Id) [lindex $ElementInfo 4]
    #set N3(Id) [lindex $ElementInfo 5]
    #set N4(Id) [lindex $ElementInfo 6]
    #set N5(Id) [lindex $ElementInfo 7]
    #set N6(Id) [lindex $ElementInfo 8]
    #set N7(Id) [lindex $ElementInfo 9]
    #set N8(Id) [lindex $ElementInfo 10]

    ## Obtaining nodes coordinates
    #set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    #set N1(x) [lindex $NCoord 0]
    #set N1(y) [lindex $NCoord 1]
    #set N1(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    #set N2(x) [lindex $NCoord 0]
    #set N2(y) [lindex $NCoord 1]
    #set N2(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    #set N3(x) [lindex $NCoord 0]
    #set N3(y) [lindex $NCoord 1]
    #set N3(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    #set N4(x) [lindex $NCoord 0]
    #set N4(y) [lindex $NCoord 1]
    #set N4(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N5(Id)] 0]
    #set N5(x) [lindex $NCoord 0]
    #set N5(y) [lindex $NCoord 1]
    #set N5(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N6(Id)] 0]
    #set N6(x) [lindex $NCoord 0]
    #set N6(y) [lindex $NCoord 1]
    #set N6(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N7(Id)] 0]
    #set N7(x) [lindex $NCoord 0]
    #set N7(y) [lindex $NCoord 1]
    #set N7(z) [lindex $NCoord 2]
    #set NCoord [lindex [GiD_Info Coordinates $N8(Id)] 0]
    #set N8(x) [lindex $NCoord 0]
    #set N8(y) [lindex $NCoord 1]
    #set N8(z) [lindex $NCoord 2]

    ## Computing element lengths
    #set lx [expr { sqrt( (0.25*($N2(x)+$N6(x)+$N3(x)+$N7(x)-$N1(x)-$N5(x)-$N4(x)-$N8(x)))**2 + (0.25*($N2(y)+$N6(y)+$N3(y)+$N7(y)-$N1(y)-$N5(y)-$N4(y)-$N8(y)))**2 + (0.25*($N2(z)+$N6(z)+$N3(z)+$N7(z)-$N1(z)-$N5(z)-$N4(z)-$N8(z)))**2 ) }]
    #set ly [expr { sqrt( (0.25*($N3(x)+$N4(x)+$N7(x)+$N8(x)-$N1(x)-$N2(x)-$N5(x)-$N6(x)))**2 + (0.25*($N3(y)+$N4(y)+$N7(y)+$N8(y)-$N1(y)-$N2(y)-$N5(y)-$N6(y)))**2 + (0.25*($N3(z)+$N4(z)+$N7(z)+$N8(z)-$N1(z)-$N2(z)-$N5(z)-$N6(z)))**2 ) }]
    #set lz [expr { sqrt( (0.25*($N5(x)+$N6(x)+$N7(x)+$N8(x)-$N1(x)-$N2(x)-$N3(x)-$N4(x)))**2 + (0.25*($N5(y)+$N6(y)+$N7(y)+$N8(y)-$N1(y)-$N2(y)-$N3(y)-$N4(y)))**2 + (0.25*($N5(z)+$N6(z)+$N7(z)+$N8(z)-$N1(z)-$N2(z)-$N3(z)-$N4(z)))**2 ) }]

    #if {$lz <= $lx} {
        #if {$lz <= $ly} {
            ## lz <= lx && lz <= ly
            #set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id) $N5(Id) $N6(Id) $N7(Id) $N8(Id)"
        #} else {
            ## ly < lz <= lx
            #set Connectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
        #}
    #} elseif {$ly <= $lx} {
        ## ly <= lx < lz
        #set Connectivities "$N1(Id) $N4(Id) $N8(Id) $N5(Id) $N2(Id) $N3(Id) $N7(Id) $N6(Id)"
    #} else {
        ## lx < lz && lx < ly
        #set Connectivities "$N1(Id) $N5(Id) $N6(Id) $N2(Id) $N4(Id) $N8(Id) $N7(Id) $N3(Id)"
    #}

    #return $Connectivities
#}

#-------------------------------------------------------------------------------

proc Line2D2Connectivities { ElemId } {

    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4]"
}

#-------------------------------------------------------------------------------

proc TriangleInterface3D4Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 4]\
    [lindex $ElementInfo 5]"

    #~ set N1(Id) [lindex $ElementInfo 3]
    #~ set N2(Id) [lindex $ElementInfo 4]
    #~ set N3(Id) $N2(Id)
    #~ set N4(Id) [lindex $ElementInfo 5]

    #~ # Obtaining nodes coordinates
    #~ set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    #~ set N1(x) [lindex $NCoord 0]
    #~ set N1(y) [lindex $NCoord 1]
    #~ set N1(z) [lindex $NCoord 2]
    #~ set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    #~ set N2(x) [lindex $NCoord 0]
    #~ set N2(y) [lindex $NCoord 1]
    #~ set N2(z) [lindex $NCoord 2]
    #~ #set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    #~ set N3(x) $N2(x)
    #~ set N3(y) $N2(y)
    #~ set N3(z) $N2(z)
    #~ set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    #~ set N4(x) [lindex $NCoord 0]
    #~ set N4(y) [lindex $NCoord 1]
    #~ set N4(z) [lindex $NCoord 2]

    #~ # Computing element lengths
    #~ set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 + (0.5*($N2(z)+$N3(z)-$N1(z)-$N4(z)))**2 ) }]
    #~ set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 + (0.5*($N3(z)+$N4(z)-$N1(z)-$N2(z)))**2 ) }]

    #~ if {$ly <= $lx} {
        #~ set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    #~ } else {
        #~ set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    #~ }

    #~ return $Connectivities
}

#-------------------------------------------------------------------------------

proc QuadrilateralInterface3D4Connectivities { ElemId } {

    # Obtaining element nodes
    set ElementInfo [GiD_Mesh get element $ElemId]
    #ElementInfo: <layer> <elemtype> <NumNodes> <N1> <N2> ...
    return "[lindex $ElementInfo 3] [lindex $ElementInfo 4] [lindex $ElementInfo 5]\
    [lindex $ElementInfo 6]"

    #~ set N1(Id) [lindex $ElementInfo 3]
    #~ set N2(Id) [lindex $ElementInfo 4]
    #~ set N3(Id) [lindex $ElementInfo 5]
    #~ set N4(Id) [lindex $ElementInfo 6]

    #~ # Obtaining nodes coordinates
    #~ set NCoord [lindex [GiD_Info Coordinates $N1(Id)] 0]
    #~ set N1(x) [lindex $NCoord 0]
    #~ set N1(y) [lindex $NCoord 1]
    #~ set N1(z) [lindex $NCoord 2]
    #~ set NCoord [lindex [GiD_Info Coordinates $N2(Id)] 0]
    #~ set N2(x) [lindex $NCoord 0]
    #~ set N2(y) [lindex $NCoord 1]
    #~ set N2(z) [lindex $NCoord 2]
    #~ set NCoord [lindex [GiD_Info Coordinates $N3(Id)] 0]
    #~ set N3(x) [lindex $NCoord 0]
    #~ set N3(y) [lindex $NCoord 1]
    #~ set N3(z) [lindex $NCoord 2]
    #~ set NCoord [lindex [GiD_Info Coordinates $N4(Id)] 0]
    #~ set N4(x) [lindex $NCoord 0]
    #~ set N4(y) [lindex $NCoord 1]
    #~ set N4(z) [lindex $NCoord 2]

    #~ # Computing element lengths
    #~ set lx [expr { sqrt( (0.5*($N2(x)+$N3(x)-$N1(x)-$N4(x)))**2 + (0.5*($N2(y)+$N3(y)-$N1(y)-$N4(y)))**2 + (0.5*($N2(z)+$N3(z)-$N1(z)-$N4(z)))**2 ) }]
    #~ set ly [expr { sqrt( (0.5*($N3(x)+$N4(x)-$N1(x)-$N2(x)))**2 + (0.5*($N3(y)+$N4(y)-$N1(y)-$N2(y)))**2 + (0.5*($N3(z)+$N4(z)-$N1(z)-$N2(z)))**2 ) }]

    #~ if {$ly <= $lx} {
        #~ set Connectivities "$N1(Id) $N2(Id) $N3(Id) $N4(Id)"
    #~ } else {
        #~ set Connectivities "$N4(Id) $N1(Id) $N2(Id) $N3(Id)"
    #~ }

    #~ return $Connectivities
}


#--------------------------------------------------------------------------------------------------------------------------------------------------------------
proc WriteCableElementSubmodelPart {FileVar CondName CableElementDict} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            #puts CableElementDict
            #set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements]
            set CableList [dict get $CableElementDict $i]
            
            #set model_part_affected 0
            #set affected_elements [lindex $CableList 4]
            
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Tables
            puts $MyFileVar "  Begin SubModelPartTables"
            puts $MyFileVar "  End SubModelPartTables"
            # Nodes 
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            puts $MyFileVar "    [lindex $CableList 2]"
            puts $MyFileVar "    [lindex $CableList 3]"
            puts $MyFileVar "  End SubModelPartNodes"    
            # Elements
            puts $MyFileVar "  Begin SubModelPartElements"
            puts $MyFileVar "    [lindex $CableList 0]"
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            puts $MyFileVar "  Begin SubModelPartConditions"
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
            
            # for {set j 0} {$j < [llength $Entities]} {incr j} {

                # if {[lindex $Entities $j] in $affected_elements} {
                    # set model_part_affected 1
                # } else {
                    # puts $MyFileVar "    [lindex $Entities $j]"
                # }        
            # }
            # if {$model_part_affected} {
                # puts $MyFileVar "[lindex $CableList 0] [lindex $CableList 2] [lindex $CableList 3]"    
            # }
            # puts $MyFileVar "  End SubModelPartElements"
        }
    }
}




proc WriteElementSubmodelPart {FileVar CondName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Tables
            puts $MyFileVar "  Begin SubModelPartTables"
            puts $MyFileVar "  End SubModelPartTables"
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"
            # Elements
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements]
            puts $MyFileVar "  Begin SubModelPartElements"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            puts $MyFileVar "  Begin SubModelPartConditions"
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
        }
    }
}

#-------------------------------------------------------------------------------

proc WritePropUnionElementSubmodelPart {FileVar PropUnionElementList} {
    upvar $FileVar MyFileVar

    puts $MyFileVar "Begin SubModelPart PropagationUnion_3d_6"
    # Tables
    puts $MyFileVar "  Begin SubModelPartTables"
    puts $MyFileVar "  End SubModelPartTables"
    # Nodes
    set Entities [GiD_EntitiesGroups get "PropagationUnion_3d_6" nodes]
    puts $MyFileVar "  Begin SubModelPartNodes"
    for {set i 0} {$i < [llength $Entities]} {incr i} {
        puts $MyFileVar "    [lindex $Entities $i]"
    }
    puts $MyFileVar "  End SubModelPartNodes"
    # Elements
    puts $MyFileVar "  Begin SubModelPartElements"
    for {set i 0} {$i < [llength $PropUnionElementList]} {incr i} {
        puts $MyFileVar "    [lindex $PropUnionElementList $i]"
    }
    puts $MyFileVar "  End SubModelPartElements"
    # Conditions
    puts $MyFileVar "  Begin SubModelPartConditions"
    puts $MyFileVar "  End SubModelPartConditions"
    puts $MyFileVar "End SubModelPart"
    puts $MyFileVar ""
}

#-------------------------------------------------------------------------------

proc WriteConstraintSubmodelPart {FileVar CondName TableDict} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Tables
            set TableList [dict get $TableDict [lindex [lindex $Groups $i] 1] TableList]
            puts $MyFileVar "  Begin SubModelPartTables"
            for {set j 0} {$j < [llength $TableList]} {incr j} {
                puts $MyFileVar "    [lindex $TableList $j]"
            }
            puts $MyFileVar "  End SubModelPartTables"
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"
            # Elements
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements]
            puts $MyFileVar "  Begin SubModelPartElements"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            puts $MyFileVar "  Begin SubModelPartConditions"
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteExcavationSubmodelPart {FileVar CondName CableElementDict} {
 set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar


        set affected_elements {}
        set added_entities {}
        dict for {k CableList} $CableElementDict {
            set affected_elements [list {*}$affected_elements {*}[lindex $CableList 4]]
        }
        
        for {set i 0} {$i < [llength $Groups]} {incr i} {

            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"

            puts $MyFileVar "  Begin SubModelPartElements"
            dict for {k CableList} $CableElementDict {
                set model_part_affected 0

                # Elements
                set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] elements]
                        
                for {set j 0} {$j < [llength $Entities]} {incr j} {
                    # check if excavation is applied on cable or truss
                    if {[lindex $Entities $j] in $affected_elements} {
                        if {[lindex $Entities $j] in [lindex $CableList 4]} {
                            set model_part_affected 1
                        }
                    # write elements only when excavation is not applied on cable or truss
                    } elseif {!([lindex $Entities $j] in $added_entities)} {
                        puts $MyFileVar "    [lindex $Entities $j]"
                        lappend added_entities  [lindex $Entities $j]
                    }
                }
                if {$model_part_affected} {
                    puts $MyFileVar "    [lindex $CableList 0]"
                }
            }
            puts $MyFileVar "  End SubModelPartElements"

            # Conditions
            puts $MyFileVar "  Begin SubModelPartConditions"
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
            
        }
    }
}


#-------------------------------------------------------------------------------

proc WriteLoadSubmodelPart {FileVar CondName TableDict ConditionDict} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Tables
            set TableList [dict get $TableDict [lindex [lindex $Groups $i] 1] TableList]
            puts $MyFileVar "  Begin SubModelPartTables"
            for {set j 0} {$j < [llength $TableList]} {incr j} {
                puts $MyFileVar "    [lindex $TableList $j]"
            }
            puts $MyFileVar "  End SubModelPartTables"
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"
            # Elements
            puts $MyFileVar "  Begin SubModelPartElements"
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            set ConditionList [dict get $ConditionDict [lindex [lindex $Groups $i] 1]]
            puts $MyFileVar "  Begin SubModelPartConditions"
            for {set j 0} {$j < [llength $ConditionList]} {incr j} {
                puts $MyFileVar "    [lindex $ConditionList $j]"
            }
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
        }
    }
}

#-------------------------------------------------------------------------------

proc WriteRecordResultSubmodelPart {FileVar CondName} {
    set Groups [GiD_Info conditions $CondName groups]
    if {[llength $Groups]>0} {
        upvar $FileVar MyFileVar

        for {set i 0} {$i < [llength $Groups]} {incr i} {
            puts $MyFileVar "Begin SubModelPart [lindex [lindex $Groups $i] 1]"
            # Tables
            puts $MyFileVar "  Begin SubModelPartTables"
            puts $MyFileVar "  End SubModelPartTables"
            # Nodes
            set Entities [GiD_EntitiesGroups get [lindex [lindex $Groups $i] 1] nodes]
            puts $MyFileVar "  Begin SubModelPartNodes"
            for {set j 0} {$j < [llength $Entities]} {incr j} {
                puts $MyFileVar "    [lindex $Entities $j]"
            }
            puts $MyFileVar "  End SubModelPartNodes"
            # Elements
            puts $MyFileVar "  Begin SubModelPartElements"
            puts $MyFileVar "  End SubModelPartElements"
            # Conditions
            puts $MyFileVar "  Begin SubModelPartConditions"
            puts $MyFileVar "  End SubModelPartConditions"
            puts $MyFileVar "End SubModelPart"
            puts $MyFileVar ""
        }
    }
}
