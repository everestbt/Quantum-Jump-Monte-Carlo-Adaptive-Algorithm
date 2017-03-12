# Quantum-Jump-Monte-Carlo-Adaptive-Algorithm
An easily modifiable quantum jump Monte Carlo code that is designed to use more memory for faster results.

It is best used to linear times, and the base algorithm does not accept custom time lists.

To run the QJMCAA  the following must be passed to the function QJMCAA:

- QJMCAA.Settings() - this contains the time, the number of trajectories, the number of points and the accuracy Magnitude (see definition below)
- QJMCAA.SavingSettings() - this controls how your data will be saved. You can give it parameters to include in the name, add expectation values to save (this is a first come first served method) and the accuracy. This is also where you set Histograms to be one or not, and again when they are added it is first come first served. Note you must add the expectation to be saved before the histogram.
- H - this is the Hamiltonian of your system in scipy.sparse.csc_matrix format.
- jumpOps - a list [] of your jump operators in scipy.sparse.csc_matrix format.
- eOps - a list [] of your expecation operators in scipy.sparse.csc_matrix format. Note unless they are added in the SavingSettings they will not be saved and their order must align.
- psi0 - your initial state in numpy.ndarray format.

Please view the twoLevel.py file for a clear example of how to use it

The accuracy magnitude (K) defines the smallest time that the algorithm will consider. It does this by taking the dt defined by the parameters T and the number of points (M) (i.e dt = T / (M-1)) and it will consider time values which are K orders of magnitude smaller. I suggest you define your accuracy magnitude as a function of T and M and the largest rate possible in the system

Dependencies:
numpy
scipy
progressbar
datetime

To run examples:
qutip
matplotlib


To run tests:
qutip
unittest
nose
