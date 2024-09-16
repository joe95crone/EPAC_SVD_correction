#!/bin/bash 

# Runs using correction.py -$4 seed -${11} trial -${13} beam dist no
elegant correct.ele -macro=prefix=$4_${11}_${13},lattice=$1,beamline=$2,energy=$3,seed=$4,InitErr=$5,InitErrXYP=$6,PMQErr=$7,PMQFSEErr=$8,EMQErr=$9,param=${10},beam_path=${12}

# correction run
sdds2plaindata track_files/correct/$4_${11}_${13}.*.sig correction_run/sig_correct_$4_${11}_${13}.dat -col=s -col=Sx -col=Sy -separator=" "
sdds2plaindata track_files/correct/$4_${11}_${13}.*.cen correction_run/cen_correct_$4_${11}_${13}.dat -col=s -col=Cx -col=Cy -col=Cxp -col=Cyp -separator=" "
sdds2plaindata track_files/correct/$4_${11}_${13}.correct.bpmcen correction_run/$4_${11}_${13}_bpm.dat -col=s -col=Cx -col=Cy -separator=" "

sdds2plaindata track_files/correct/$4_${11}_${13}.correct.bpmcen bpm_files/correct/$4_${11}_${13}_bpm.dat -col=s -col=Cx -col=Cy -separator=" "
