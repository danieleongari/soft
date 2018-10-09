Split-Operator Fourier Transform (SOFT) program
for the course [Computational Methods in Molecular Quantum Mechanics CH-452](http://edu.epfl.ch/coursebook/en/computational-methods-in-molecular-quantum-mechanics-CH-452).

Original program written by Sara Bonella.

### HOW TO PLAY WITH IT:
- fork it
- clone it from your repository
- modify SOFT.c, e.g.,
     - grid and time step settings
     - initial position and speed
     - mass
     - potential type (barrier, morse, harmonic)
- compile it with $make
- run SOFT.x
- visualize the output .dat using gnuplot (e.g., $gnuplot psi.gnuplot)
- report any bug!

### EXAMPLE:
This is a gaussian wave packet, with a starting position and momentum, hitting a rectangular potential barrier.

![](README.gif)
