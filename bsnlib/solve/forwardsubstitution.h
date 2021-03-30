
#ifndef FORWARDSUBSTITUTION_H
#define FORWARDSUBSTITUTION_H

namespace bsnlib {

namespace solve {

/** Forward substitution */
template <typename S> struct Solve<S, SolveMethod::FORWARD> {
    template <typename A, typename B, typename C> static void solve(const A& M, B& x, const C& y) {
        const int R = internal::traits<A>::rows;
        S tmp = static_cast<S>(0);

        for (int i=0; i<R; ++i){
            for (int j=0; j<i; ++j){
                tmp += M(i, j) * x(j, 0);
            }

            x(i, 0) = (y(i, 0) - tmp) / M(i, i);
            tmp = static_cast<S>(0);
        }
    }
};


} // solve

} // bsnlib

#endif // FORWARDSUBSTITUTION_H
