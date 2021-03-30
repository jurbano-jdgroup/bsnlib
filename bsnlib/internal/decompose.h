
#ifndef DECOMPOSE_H
#define DECOMPOSE_H

namespace bsnlib {

namespace solve {

/** decompose class definition */
template <typename S, int M, int R, int C> class Decompose;

/** decompose method definition */
struct DecompositionMethod {
    enum {
        LU = 0,
        PLU = 1,
        LLT = 2,
        LDLT = 3,
        QR = 4,
        SVD = 5
    };
};

} // solve

} // bsnlib

#endif // DECOMPOSE_H
