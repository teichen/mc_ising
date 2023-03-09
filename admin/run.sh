#!/bin/bash 
#PBS -l walltime=999:00:00
#PBS -q default
#PBS -N mcapbs
#PBS -l nodes=1:ppn=24

SCR_DIR=/tmp/scr
HOME_DIR=/home

mkdir $SCR_DIR

cd $SCR_DIR

rest=0
logging=0
lambda=1
tsteps=100

time OMP_NUM_THREADS=4 OMP_STACKSIZE=100M OMP_DYNAMIC=TRUE ./glauber_ising $rest $logging $lambda $tsteps > log.txt &
cp ./n.dat ./log.txt $HOME_DIR

cd $HOME_DIR
rm -rf $SCR_DIR

