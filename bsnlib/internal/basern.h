
#ifndef BASERN_H
#define BASERN_H

//#include "bsnlib/internal/traits.h"
#include "assert.h"

namespace bsnlib {

namespace internal {

/** Base Rn definition */
template <typename S, int N, int Idx> class RnBase;

/** Base Rn class */
template <typename S, int N, int Idx> struct traits<RnBase<S, N, Idx>> {
    enum{
        rows = N,
        cols = 1
    };
};

/** Base Rn implementation */
template <typename S, int N, int Idx> class RnBase {
public:

    /** default constructor */
    RnBase () {
        static_assert (Idx<=N, "Rn base element, idx grather than N");
    }

    /** entries selection operator */
    S operator()(unsigned int i, unsigned int j) const {
        assert(i<N && j==0);
        return i==(Idx-1) ? static_cast<S>(1) : static_cast<S>(0);
    }
};

} // internal

} // bsnlib

#endif // BASERN_H
