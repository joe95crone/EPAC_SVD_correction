#!/bin/bash

rm track_files/temp/*

elegant lat_gen.ele -macro=lattice=$1,beamline=$2,energy=$3

# getting a baseline for each error run
sdds2plaindata track_files/temp/*.twi baseline/twiss_base.dat -col=s -col=beta? -col=alpha? -col=psi? -col=eta? -separator=" "
sdds2plaindata track_files/temp/*.sig baseline/sig_base.dat -col=s -col=Sx -col=Sy -col=beta?Beam -col=alpha?Beam -separator=" "
sdds2plaindata track_files/temp/*.cen baseline/cen_base.dat -col=s -col=Cx -col=Cy -col=Cxp -col=Cyp -separator=" "
sdds2plaindata track_files/temp/*.bpmcen baseline/bpm_base.dat -col=s -col=Cx -col=Cy -separator=" "
