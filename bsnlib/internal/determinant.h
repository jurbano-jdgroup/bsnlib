
#ifndef DETERMINANT_H
#define DETERMINANT_H

//#include "bsnlib/internal/traits.h"

namespace bsnlib {

namespace internal {

/** Determinant implementation */
template <typename S, int R> struct Det {};


/** 2x2 matrix determinant  */
template <typename S> struct Det<S, 2> {
    template <typename M> static S calculate(const M& m) {
        return m(0, 0)*M(1, 1) - M(0, 1)*M(1, 0);
    }
};

/** 3x3 matrix determinant */
template <typename S> struct Det<S, 3> {
    template <typename M> static S calculate(const M& m) {
        const S x_in = m(1, 1)*m(2, 2) - m(1, 2)*m(2, 1);
        const S y_in = m(1, 0)*m(2, 2) - m(1, 2)*m(2, 0);
        const S z_in = m(1, 0)*m(2, 1) - m(1, 1)*m(2, 0);

        return m(0,0)*x_in - m(0, 1)*y_in + m(0,2)*z_in;
    }
};

/** Determinant struct */
template <typename S> struct Determinant {
    template <typename M> static S calculate(const M& m) {
        /* check the matrix dimensions */
        static_assert (traits<M>::rows<=3,"Matrix dimension greather than 3 is not implemented");
        /* check the S type */
        static_assert(traits<S>::rows==1 && traits<S>::cols==1, "Scalar type dimension error");
        /* check the matrix square */
        static_assert(int(traits<M>::rows)==int(traits<M>::cols), "Determinant on non square matrix");
        return Det<S, traits<M>::rows>::template calculate<M>(m);
    }
};

} // internal

} // bsnlib

#endif // DETERMINANT_H
