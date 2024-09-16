#!/bin/bash

# remove previous lattice
rm temp.lte
rm track_files/temp/*

# remove response matrix
rm bpm_files/response/*
rm track_files/response/*
rm response_matrices/*

# remove error runs
rm bpm_files/error/*
rm track_files/error/*
rm track_files/error/logs/*
rm error_run/*

# remove correction runs
rm bpm_files/correct/*
rm track_files/correct/*
rm correction_run/*
 
