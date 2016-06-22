# mdsim
Molecular Dynamics: apply the linked cell algorithm to our simulation
of two colliding blocks consisting of several thousan particles

In order to run this programm use the following command:
```bash
$ ./mdsim [parameter file] [data file]
```

If parameters are chosen wrong, velocities or forces can get too big -> Segmentation fault

Possible parameter files:
blocks.par
mini.par

Possible data files:
blocks-large.dat
blocks-medium.dat
blocks-small.dat
mini.dat

This file includes:
  - Makefile
  - mdsim.cpp
  - README
