#ifndef _HELPFUL_GLOBAL_MACROS_H
#define _HELPFUL_GLOBAL_MACROS_H

// This file just contains nice macros 
// that can be used everywhere without 
// breaking anything.

// Branchless ternary operator to be used in lieu of
// the branching code: cond_a ? a : b
#define _BRANCHLESS_TERNARY( cond_a, a, b ) ( (cond_a) * (a) + !(cond_a) * (b) )


#endif /* _HELPFUL_GLOBAL_MACROS_H */