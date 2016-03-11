#!/bin/bash
#MSUB -l nodes=1:ppn=1
#MSUB -l walltime=01:00:00:00
#MSUB -N 1009889_standard_[]
#MSUB -A p20254
#MSUB -o jobout.txt
#MSUB -e joberr.txt

module load mpi/openmpi-1.6.3-intel2011.3
module load intel/2011.3
ulimit -s unlimited

cd $PBS_O_WORKDIR
NPROCS=`wc -l < $PBS_NODEFILE`
#running on quest

gunzip -f CHGCAR.gz &> /dev/null
date +%s
ulimit -s unlimited

 /projects/b1004/bin/vasp_53_serial  > stdout.txt 2> stderr.txt

gzip -f CHGCAR OUTCAR PROCAR ELFCAR
rm -f WAVECAR vasprun.xml
date +%s
