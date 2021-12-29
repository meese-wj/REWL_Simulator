#ifndef SQUARE_2D_NEAREST_NEIGHBOR_TESTER_H
#define SQUARE_2D_NEAREST_NEIGHBOR_TESTER_H

// This is a class that tests whether the 
// address book for the nearest neighbors works
// on a 2D square lattice.

// This test fills the neighbors in the order of: (x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)

#include "../square_2d_nearest_neighbors.hpp"
#include <vector>
#include <iostream>

class square_2d_nearest_neighbor_tester 
{
public:
    using neighbor_answers = std::vector<std::array<std::uint32_t, 4> >;
    square_2d_nearest_neighbor_tester(/* args */);
    int run_test() const;
    ~square_2d_nearest_neighbor_tester();
private:
    neighbor_answers three_by_three = { 
                                        {2, 1, 6, 3},    // 0 = (0, 0)
                                        {0, 2, 7, 4},    // 1 = (1, 0) 
                                        {1, 0, 8, 5},    // 2 = (2, 0) 
                                        {5, 4, 0, 6},    // 3 = (0, 1) 
                                        {3, 5, 1, 7},    // 4 = (1, 1) 
                                        {4, 3, 2, 8},    // 5 = (2, 1) 
                                        {8, 7, 3, 0},    // 6 = (0, 2) 
                                        {6, 8, 4, 1},    // 7 = (1, 2) 
                                        {7, 6, 5, 2},    // 8 = (2, 2) 
                                                        };
    neighbor_answers three_by_two = {
                                        {2, 1, 3, 3},    // 0 = (0, 0)
                                        {0, 2, 4, 4},    // 1 = (1, 0) 
                                        {1, 0, 5, 5},    // 2 = (2, 0) 
                                        {5, 4, 0, 0},    // 3 = (0, 1) 
                                        {3, 5, 1, 1},    // 4 = (1, 1) 
                                        {4, 3, 2, 2},    // 5 = (2, 1)
                                                        };
    neighbor_answers two_by_three = {
                                        {1, 1, 4, 2},    // 0 = (0, 0)
                                        {0, 0, 5, 3},    // 1 = (1, 0) 
                                        {3, 3, 0, 4},    // 2 = (0, 1) 
                                        {2, 2, 1, 5},    // 3 = (1, 1) 
                                        {5, 5, 2, 0},    // 4 = (2, 0) 
                                        {4, 4, 3, 1},    // 5 = (2, 1)
                                                        };
};

square_2d_nearest_neighbor_tester::square_2d_nearest_neighbor_tester(/* args */)
{
}

square_2d_nearest_neighbor_tester::~square_2d_nearest_neighbor_tester()
{
}

int square_2d_nearest_neighbor_tester::run_test() const
{
    int test_passed = true;

    Square_2D_Nearest_Neighbors book1( 3, 3 );
    book1.initialize();
    for ( std::uint32_t site = 0; site != 3 * 3; ++site )
        for (std::uint32_t nn = 0; nn != 4; ++nn )
            test_passed *= ( three_by_three[site][nn] == book1.neighbor(site, nn) );

    Square_2D_Nearest_Neighbors book2( 3, 2 );
    book2.initialize();
    for ( std::uint32_t site = 0; site != 3 * 2; ++site )
        for (std::uint32_t nn = 0; nn != 4; ++nn )
            test_passed *= ( three_by_two[site][nn] == book2.neighbor(site, nn) );

    Square_2D_Nearest_Neighbors book3( 2, 3 );
    book3.initialize();
    for ( std::uint32_t site = 0; site != 2 * 3; ++site )
        for (std::uint32_t nn = 0; nn != 4; ++nn )
            test_passed *= ( two_by_three[site][nn] == book3.neighbor(site, nn) );

    return test_passed;
}



#endif /* SQUARE_2D_NEAREST_NEIGHBOR_TESTER_H */