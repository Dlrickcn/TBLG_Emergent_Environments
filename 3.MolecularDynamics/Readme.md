# Phonon total density of state curves calculated from molecular dynamics simulations for 2.13º rotated TBLG 

This study shows the usage of fix-phonon to calculate the dynamical matrix as well as phonon dispersion curve for TBLG based on a Tersoff potential.

The files under this directory:
1) Graphene.bin.20000000  : last output binary file by fix-phonon
2) data.pos               : LAMMPS input file
3) pdos.dat              : phonon DOS  data from Graphene.bin.20000000
4) in.graphene            : LAMMPS input file
5) map.in                 : LAMMPS input file for fix-phonon
6) pdos.png               : figure of phonon DOS curves
7) plot.pdos              : gnuplot script to generate pdos.eps
8) pdos.gnuplot           : gnuplot script to generate pdos.eps (auto generated)

To run this example, simply invoke: 
-> lmp -in in.graphene -screen none

Once done, one can use the auxiliary analysing code "phana" to obtain "pdisp.dat"
-> phana Graphene.bin.20000000 < in.dos

And then use the gnuplot script file "plot.dos" to generate pdisp.eps:
-> gnuplot plot.pdos

The resultant ``pdos.png'' shows the measured phonon dispersion.

**Our work is used the details of the code in the [link](https://github.com/lingtikong/fix-phonon/tree/master)**


