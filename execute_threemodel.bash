#!/bin/bash

#Compile:
ifort re_const.f90 -o re_const 
ifort re_gau.f90 -o re_gau
ifort re_pl.f90 -o re_pl

#Execute:
nohup ./re_const > const_relog.txt 2>&1 &
nohup ./re_gau > gau_relog.txt 2>&1 &
nohup ./re_pl > pl_relog.txt 2>&1 &

