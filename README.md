qualifier
=========

My PhD take home qualifying exam.  2 Computational projects.

=========
Basin Hopping:

Main is in basinhop.cpp, as is the Monte Carlo step.  The total energy calculation and minimization methods (A&B) are in the localMin.cpp file.  The structure 'state' stores all relevant information for a cluster including:
-Atomic Coordinates
-Center of mass
-Energy
-Storage space for calcuations

Powell's method was taken from Numerical Recipes this includes the files:
f1dim.c
linmin.c
nr.h
nrutil.c
nrutil.h
powell.c
mnbrak.c

Several python scripts were written for visualization the results and making graphs to answer posted questions.  The script viz.py takes as its only arguement the location of a log file, e.g.:
./viz.py ../basinhopping_data/finalstate_N104_A_0.dat.
Which prints total energy and plots the lowest energy atomic structure with bonded neighbors highlighted by a line.
The others can be run as follows
./viz12.py <# atoms>
./viz3.py <#atoms>
./viz4.py <#atoms>
./viz5.py -no args-qualifier
=========