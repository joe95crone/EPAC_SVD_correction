#!/bin/bash

rm track_files/simple/*
rm watch_files/simple/*

elegant simple_run.ele

# Matrix Twiss 
#sddsplot track_files/simple/*.twi -layout=2,1 \
# -col=s,beta? -graphic=line,type=1,vary -newPanel -title=beta -endPanel -legend \
#-col=s,alpha? -graphic=line,type=3,vary -newPanel -title=alpha -endPanel -topline=matrix -legend 

# Matrix Twiss + Dispersion 
sddsplot track_files/simple/*.twi -TopTitle -layout=3,1 \
 -col=s,beta? -graphic=line,type=1,vary -newPanel -title=beta -endPanel -legend\
 -col=s,alpha? -graphic=line,type=3,vary -newPanel -title=alpha -endPanel -legend\
 -col=s,eta? -graphic=line,type=5,vary -newPanel -title=dispersion -endPanel -legend


# Tracked Twiss 
#sddsplot track_files/simple/*.sig -layout=2,1 \
# -col=s,beta?Beam -graphic=line,type=1,vary -newPanel -title=beta -endPanel \
#-col=s,alpha?Beam -graphic=line,type=3,vary -newPanel -title=alpha -endPanel -topline=tracked 

# FBPIC Beam (Nicolas Bourgeois)
#sdds2plaindata track_files/*.sig Simple_Analysis/FFAGFULL_NB_cut_scaled_sig.dat -col=s -col=s16 -col=s6 -col=Sx -col=Sy -col=beta?Beam -col=alpha?Beam -separator=" "
#sdds2plaindata track_files/*.out Simple_Analysis/FFAGFULL_NB_cut_scaled_out.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" " 
#sdds2plaindata FBPIC_Beam/Saved_Beams/*scale.sdds Simple_Analysis/FFAGFULL_NB_cut_scaled_in.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" "

sdds2plaindata track_files/simple/*.twi simple_run/twiss_matrix.dat -col=s -col=beta? -col=alpha? -col=psi? -col=eta? -separator=" "
sdds2plaindata track_files/simple/*.sig simple_run/twiss_tracked.dat -col=s -col=Sx -col=Sy -col=beta?Beam -col=alpha?Beam -separator=" "
sdds2plaindata track_files/simple/*.out simple_run/beam_tracked.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" "

#sdds2plaindata track_files/simple/*.bun simple_run/Pseudo_CUT.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" "
#sdds2plaindata ../FBPIC_Beam/Saved_Beams/OF2_EPAC_KDE_1GeV.sdds simple_run/FBPIC_FULL.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" "

# WATCH POINTS THROUGHOUT LATTICE (11 for 4 PMQ 3 EMQ)
#for i in {1..11}
#do
#	if [ $i -lt 9 ]
#	then
#		sdds2plaindata watch_files/simple/*00$i.w1 simple_run/WATCH/beam_tracked_$i.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" "
#	else
#		sdds2plaindata watch_files/simple/*$i.w1 simple_run/WATCH/beam_tracked_$i.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" "
#	fi
#done

# WATCH POINT SLICES IN SPECTROMETER
# Logic required to distinguish between *11.wspec and *01.wspec
#for i in {1..15}
#do
#	if [ $i -lt 9 ]; 
#	then 
#		sdds2plaindata watch_files/simple/*00$i.wspec simple_run/WATCH_SPEC/beam_tracked_$i.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" ";
#	else
#		sdds2plaindata watch_files/simple/*$i.wspec simple_run/WATCH_SPEC/beam_tracked_$i.dat -col=x -col=xp -col=y -col=yp -col=t -col=p -separator=" ";
#	fi
#done
