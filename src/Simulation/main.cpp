#include <stdio.h>

#include "rewl_simulation.hpp"

int main(const int argc, const char * argv[])
{
    if (argc > 1)
    {
        printf("\nThis simulation %s takes no command line arguments.", argv[0]);
        printf("\nReturning error value.\n\n");
        return 1;
    }


    REWL_simulation * simulation = new REWL_simulation();
    

    delete simulation;

    return 0;
}
