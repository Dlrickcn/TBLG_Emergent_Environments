# Phonon total density of state curves calculated from molecular dynamics simulations for 2.13º rotated TBLG 

This study shows the usage of fix-phonon to calculate the dynamical matrix as well as phonon dispersion curve for TBLG based on a Tersoff potential.

The files under this directory:
 1) Graphene.bin.20000000   : last output binary file by fix-phonon
 2) data.pos               : LAMMPS input file
 3) pdisp.dat              : phonon dispersion data from Graphene.bin.20000000
 4) in.disp                : input file to get disp.dat by phana
 5) in.graphene            : LAMMPS input file
 6) log.lammps             : LAMMPS log file
 7) map.in                 : LAMMPS input file for fix-phonon
 8) pdisp.png              : figure of phonon dispersion curves
10) plot.disp              : gnuplot script to generate pdisp.eps
11) pdisp.gnuplot          : gnuplot script to generate pdisp.eps (auto generated)

To run this example, simply invoke: 
-> lmp -in in.graphene -screen none

Once done, one can use the auxiliary analysing code "phana" to obtain "pdisp.dat"
-> phana Graphene.bin.20000000 < in.disp

And then use the gnuplot script file "plot.disp" to generate pdisp.eps:
-> gnuplot plot.pdisp

The resultant ``pdisp.png'' shows the measured phonon dispersion.

**Our work is used the details of the code in the [link](https://github.com/lingtikong/fix-phonon/tree/master)**


