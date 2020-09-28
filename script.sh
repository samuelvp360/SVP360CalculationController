#!/bin/bash

#$ -pe mpi 32 #no. of processors going to using, can assign up to 128
#$ -N LAM1 #job name
#$ -cwd #change to working directory
#$ -e output2.err #specify file name for standard error
#$ -o output2.out #specify file name for standard output
# setting up Guassian environment
export g16root=/home/samuelvip
# $HOME/tmp folder already been changed by chmod -R 777 $HOME/tmp
export GAUSS_SCRDIR=/home/samuelvip/g16/scratch
source $g16root/g16/bsd/g16.profile
# export GAUSS_EXEDIR=/home/samuelvip/g16

$g16root/g16/g16 < Opt_TCE.com > Opt_TCE.log
