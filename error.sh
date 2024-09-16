#!/bin/bash

# -$10 is seed, -$11 is no. trials, -$13 is beam dist no 
elegant error.ele -macro=prefix=${10}_${11}_${13},lattice=$1,beamline=$2,energy=$3,seed=$4,InitErr=$5,InitErrXYP=$6,PMQErr=$7,PMQFSEErr=$8,EMQErr=$9,beam_path=${12} 

# Error Run
sdds2plaindata track_files/error/${10}_${11}_${13}.*.sig error_run/sig_error_${10}_${11}_${13}.dat -col=s -col=Sx -col=Sy -separator=" "
sdds2plaindata track_files/error/${10}_${11}_${13}.*.cen error_run/cen_error_${10}_${11}_${13}.dat -col=s -col=Cx -col=Cy -col=Cxp -col=Cyp -separator=" "
sdds2plaindata track_files/error/logs/${10}_${11}_${13}*.erl error_run/log_error_${10}_${11}_${13}.dat -col=ElementName -col=ElementParameter -col=ParameterValue -separator=" " 
sdds2plaindata track_files/error/${10}_${11}_${13}.*.bpmcen error_run/${10}_${11}_${13}_bpm.dat -col=s -col=Cx -col=Cy -Separator=" "

sdds2plaindata track_files/error/${10}_${11}_${13}.*.bpmcen bpm_files/error/${10}_${11}_${13}_bpm.dat -col=s -col=Cx -col=Cy -Separator=" "
