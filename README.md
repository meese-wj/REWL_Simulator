# Replica-Exchange Wang Landau Simulator
### W. Joe Meese ~ meese022@umn.edu

## Summary 
This code base is a _Replica-Exchange Wang Landau_ (REWL) simulator written to generalize the REWL method to various Hamiltonians smoothly. Once written, it will allow the user to only specify the details of a templated Hamiltonian and observables and remove the user from the nitty-gritty details of replica communication via MPI. In principle, this will allow quick development and study of different classical thermodynamic systems without the need for regular development of the
actual algorithm itself.

## Code Basics
The goal for this code base is for it to be written in `C++` and built with `cmake` across different platforms. Currently it has only been tested on Linux systems. This section will be updated as more of the code is developed.

### Code Installation
First, clone this repo directly. Suppose the relative path to the clone is then `repo_path`. Then `cd` into that directory as `cd repo_path`. Next create a build directory within the cloned source. _It is not recommended to build in-source!_ From the Unix terminal, execute the following
```bash 
# After 'cd repo_path'
mkdir build
cd build
cmake ..
make -j 
```
`cmake` will fail if an MPI implementation is not found. If it is, `cmake` will determine the correct MPI compiler and then the corresponding runtime flags. These will be outputted at the configure step of `cmake` and will be stored in the `build/CMakeCache.txt` file if the configure step is successful. 

When `make -j` is called for parallel compilation, the simulation runtime executable will be stored in the build directory as `bin/REWL_Simulator`. The runtime command to execute the code will come down to
```bash
MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_MAX_NUMPROCS MPIEXEC_PREFLAGS ./bin/REWL_Simulator MPIEXEC_POSTFLAGS <simulator-flags>
```
Currently, there are no runtime flags to be passed as `<simulator-flags>` to the simulation.

Some `cmake` options to be aware of are as follows
* `MPI_ON`
    * This engages parallelization and replica exchange. 
    * This option is deprecated and will soon be mandatory.
* `ISING2D`
    * This flag sets the Ising model on a 2d square periodic grid.
    * It is the default Hamiltonian used.
* `COLLECT_TIMINGS`
    * This flag compiles the timing functionality through the `C++` Standard Library (`std::chrono`).
* `PRINT_HISTOGRAM`
    * When set, the histograms are intermittently written to a file.
    * Currently this is _not_ thread safe and will be corrected as the parallelization is implemented.
