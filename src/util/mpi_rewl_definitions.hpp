#ifndef MPI_REWL_DEFINITIONS
#define MPI_REWL_DEFINITIONS
/* Use this file to store global
 * MPI definitions */

constexpr int REWL_MASTER_PROC = 0;

enum REWL_Tags
{
    final_energy_tag, final_logdos_tag, final_obs_tag
};

#endif
