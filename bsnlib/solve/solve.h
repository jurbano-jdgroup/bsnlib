
#ifndef SOLVE_H
#define SOLVE_H

namespace bsnlib {

namespace solve {

/**
 * Solve template struct
 * definition
 *
 */
template <typename S, int M> struct Solve;

/**
 * Solve methods
 * definition
 *
 */
struct SolveMethod{
    enum{
        BACK = 0,
        FORWARD = 1,
        DIRECT = 2
    };
};



} // solve

} // bsnlib

#endif // SOLVE_H
