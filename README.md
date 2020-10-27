# Replica-Exchange Wang Landau Simulator
### W. Joe Meese ~ meese022@umn.edu

## Summary 
This code base is a _Replica-Exchange Wang Landau_ (REWL) simulator written to generalize the REWL method to various Hamiltonians smoothly. Once written, it will allow the user to only specify the details of a templated Hamiltonian and observables and remove the user from the nitty-gritty details of replica communication via MPI. In principle, this will allow quick development and study of different classical thermodynamic systems without the need for regular development of the
actual algorithm itself.

## Code Basics
The goal for this code base is for it to be written in `C++` and built with `cmake` across different platforms. Currently it has only been tested on Linux systems. This section will be updated as more of the code is developed.

## System Requirements
The code base requires the `std::filesystem` which was officially adopted in the `C++17` standard. It is recommended that a recent enough compiler is used, or else pieces of the `filesystem` will be scattered through the standard library in `iomanip` and `std::experimental::filesystem`, and so the compilation with `make` will fail. On Intel 64-bit systems, successful compilation is achieved with `GNU 9.x` and `OpenMPI 4.x` or newer.

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
`cmake` will fail if an MPI implementation is not found. If it is, `cmake` will determine the correct MPI compiler and then the corresponding runtime flags. These will be outputted at the configure step of `cmake` and will be stored in the `build/CMakeCache.txt` file if the configure step is successful. 

When `make -j` is called for parallel compilation, the simulation runtime executable will be stored in the build directory as `bin/REWL_Simulator`. The runtime command to execute the code will come down to
```bash
MPIEXEC MPIEXEC_NUMPROC_FLAG MPIEXEC_MAX_NUMPROCS MPIEXEC_PREFLAGS ./bin/REWL_Simulator MPIEXEC_POSTFLAGS <simulator-flags>
```
Currently, there are no runtime flags to be passed as `<simulator-flags>` to the simulation.

### Build flags
Some `cmake` options to be aware of are as follows:

#### Global simulation control flags
* `MPI_ON`
    * This engages parallelization and replica exchange. 
    * This option is deprecated and will soon be mandatory.
* `COLLECT_TIMINGS`
    * This flag compiles the timing functionality through the `C++` Standard Library (`std::chrono`).
* `PRINT_HISTOGRAM`
    * When set, the histograms are intermittently written to a file.
    * Currently this is _not_ thread safe and will be corrected as the parallelization is implemented.
* `REDUCE_LOGDOS`
    * This truncates the logDoS at the same time as the energy histogram such that the ground state value of the logDoS is zero. Since only differences in the logDoS affect the Monte Carlo moves, this is allowed.
    * The idea here is to reduce the degree of numerical error saturation in calculating acceptance probabilities. When the logDoS gets too large, it will overflow the mantissa, especially when the incrementer is much smaller than the magnitude of the logDoS.
* `SAMPLE_AFTER`
    * This waits until one full round of REWL (or simply WL) has been completed before any sampling of the non-energetic observables is performed. By doing so, this guarantees that _Detailed-Balance_ is more closely approximated for these statistical observables.
    * For the second round of REWL moves with observable sampling, then the logDoS incrementer is reset back to one, and the simulation proceeds as it did before, only now it takes measurements concurrently.
* `TRAPEZOIDAL_RULE`
    * Evaluate the thermodynamic integrals with the trapezoidal rule. 
    * This will likely be absorbed into the main code in the future.
* `INDEPENDENT_WALKERS`
    * If engaged, there will be no replica exchange, nor will there be any intra-window walker averaging.
* `DIFFERENT_SEEDS`
    * Use different seeds for each walkers based on high-speed clock. If not engaged, then all seeds are identical.

#### Supported Hamiltonians
If multiple Hamiltonians are selected at build time, then `cmake` will throw a `FATAL_ERROR`.
* `ISING2D`
    * This flag sets the Ising model on a 2d square periodic grid.
    * It is the default Hamiltonian used.
* `ASHKIN_TELLER2D`
    * This flag sets the Ashkin-Teller model on a 2d square periodic grid.
    * Right now, the Ashkin-Teller model is supported within the ferromagnetic Baxter regime. 
