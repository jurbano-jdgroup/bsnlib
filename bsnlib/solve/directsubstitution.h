
#ifndef DIRECTSUBSTITUTION_H
#define DIRECTSUBSTITUTION_H

namespace bsnlib {

namespace solve {

/** Direct substitution */
template <typename S> struct Solve<S, SolveMethod::DIRECT> {
    template <typename A, typename B, typename C> static void solve(const A& M, B& x, const C& y) {
        const int R = internal::traits<A>::rows;
        for (int i=0; i<R; ++i) {

            x(i, 0) = y(i, 0) / M(i, i);
        }
    }
};


} // solve

} // bsnlib

#endif // DIRECTSUBSTITUTION_H
