#include <stdio.h>

#include <model_hamiltonians.hpp>

int main(const int argc, const char * argv[])
{
    if (argc > 1)
    {
        printf("\nThis simulation %s takes no command line arguments.", argv[0]);
        printf("\nReturning error value.\n\n");
        return 1;
    }
    return 0;
}
