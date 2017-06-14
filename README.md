# Quantum-Jump-Monte-Carlo-Adaptive-Algorithm

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/3cb553bc3e70437cb09aea28e7f98338)](https://www.codacy.com/app/everest.bt/Quantum-Jump-Monte-Carlo-Adaptive-Algorithm?utm_source=github.com&utm_medium=referral&utm_content=everestbt/Quantum-Jump-Monte-Carlo-Adaptive-Algorithm&utm_campaign=badger)

An easily modifiable quantum jump Monte Carlo code that is designed to use more memory for faster results.

It is best used to linear times, and the base algorithm does not accept custom time lists.

To run the QJMCAA the following must be passed to the function QJMCAA.QJMCRun(QJMCAA.Settings, QJMCAA.SavingSettings, H, jumpOps, eOps, psi0):

- QJMCAA.Settings() - this contains the time, the number of trajectories, the number of points and the accuracy Magnitude (see definition below)
- QJMCAA.SavingSettings() - this controls how your data will be saved. You can give it parameters to include in the name, add expectation values to save (this is a first come first served method) and the accuracy. This is also where you set Histograms to be one or not, and again when they are added it is first come first served. Note you must add the expectation to be saved before the histogram.
- H - this is the Hamiltonian of your system in scipy.sparse.csc_matrix format.
- jumpOps - a list [] of your jump operators in scipy.sparse.csc_matrix format.
- eOps - a list [] of your expecation operators in scipy.sparse.csc_matrix format. Note unless they are added in the SavingSettings they will not be saved and their order must align.
- psi0 - your initial state in numpy.ndarray format.

Please view the twoLevel.py file for a clear example of how to use it

The accuracy magnitude (K) defines the smallest time that the algorithm will consider. It does this by taking the dt defined by the parameters T and the number of points (M) (i.e dt = T / (M-1)) and it will consider time values which are K orders of magnitude smaller. I suggest you define your accuracy magnitude as a function of T, M, and the largest rate possible in the system.

The low memory method should only be used if you require the lower memory (surprise surprise) or if you want to do custom time lists (such as logarithmic time). This works by only producing the exponential matrices when they are needed, as such there should only be 1 at a time. This greatly reduces the memory but makes it much slower as this is the most costly action in the simulation. It is NOT optimised as it is not the goal of the package. If there is a greater interest in it, there are optimisations that can be done, I simply included it such that custom time lists could be done.

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

Please raise any issues you find on the github page:https://github.com/everestbt/Quantum-Jump-Monte-Carlo-Adaptive-Algorithm

If you wish to contribute please fork the code and place up there for review.

Warm regards,

Ben Everest
