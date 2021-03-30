
#ifndef BACKSUBSTITUTION_H
#define BACKSUBSTITUTION_H

namespace bsnlib {

namespace solve {

/** Backsubstitution */
template <typename S> struct Solve<S, SolveMethod::BACK> {
    template <typename A, typename B, typename C> static void solve(const A& M, B& x, const C& y) {
        const int R = internal::traits<A>::rows;

        S tmp = static_cast<S>(0);
        for (int i=0; i<R; ++i) {

            for (int j=0; j<i; ++j) {
                tmp += M(R - i - 1, R - j - 1) * x(R - j - 1, 0);
            }

            x(R - i - 1, 0) = (y(R - i - 1, 0) - tmp) / M(R - i - 1, R - i - 1);

            tmp = static_cast<S>(0);
        }
    }
};


} // solve

} // bsnlib

#endif // BACKSUBSTITUTION_H
