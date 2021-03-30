
#ifndef DOT_H
#define DOT_H

//#include "bsnlib/internal/traits.h"
//#include "bsnlib/internal/complex.h"

namespace bsnlib {

namespace internal {

/** Dot computation */
template <typename S, typename L, typename R, bool scalar_type> struct dot_computation {};

/** non scalar type dot computation */
template <typename S, typename L, typename R> struct dot_computation<S, L, R, false> {
    static S compute(const L& l, const R& r) {
        S tmp = static_cast<S>(0);
        for(unsigned int i=0; i<int(traits<L>::rows); ++i){
            tmp += l(i,0)*r(i,0);
        }
        return tmp;
    }
};

/** scalar type dot computation */
template <typename S, typename L, typename R> struct dot_computation<S, L, R, true> {
    static S compute(const L& l, const R& r) {
        S tmp = l*r;
        return tmp;
    }
};

/** Dot product implementation */
template <typename S> struct Dot {
    template <typename L, typename R> static S compute(const L& l, const R& r) {
        /* check the s type */
        static_assert(is_scalar<S>::value, "Scalar dimension error");
        /* check the r dimensions */
        static_assert(int(traits<L>::rows)==int(traits<R>::rows), "Dot l and r dimension mismatch");
        return dot_computation<S, L, R, (is_scalar<L>::value && is_scalar<R>::value)>::compute(l, r);
    }
};

} // internal

} // bsnlib

#endif // DOT_H
