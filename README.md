# Replica-Exchange Wang Landau Simulator
### W. Joe Meese ~ meese022@umn.edu

## Summary 
This code base is a _Replica-Exchange Wang Landau_ (REWL) simulator written to generalize the REWL method to various Hamiltonians smoothly. Once written, it will allow the user to only specify the details of a templated Hamiltonian and observables and remove the user from the nitty-gritty details of replica communication via MPI. In principle, this will allow quick development and study of different classical thermodynamic systems without the need for regular development of the
actual algorithm itself.

## Code Basics
This code base is written in `C++` and uses `cmake` to build across different platforms. This section will be updated as more of the code is developed.
