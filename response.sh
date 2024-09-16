#!/bin/bash 

# Runs using correction.py -$4 seed -$11 corrector no. (1 plate = 1 corrector) -${12} X or Y flag -${14} beam distribution number
elegant response.ele -macro=prefix=$4_${11}_${12}_${14},lattice=$1,beamline=$2,energy=$3,seed=$4,InitErr=$5,InitErrXYP=$6,PMQErr=$7,PMQFSEErr=$8,EMQErr=$9,param=${10},beam_path=${13}

# reads to file the BPM centroids in x, y 
sdds2plaindata track_files/response/$4_${11}_${12}_${14}.response.bpmcen bpm_files/response/$4_${11}_${12}_${14}_bpm.dat -col=s -col=Cx -col=Cy -separator=" "
