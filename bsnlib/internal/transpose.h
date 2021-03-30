#ifndef TRANSPOSE_H
#define TRANSPOSE_H

//#include <assert.h>
//#include "bsnlib/internal/traits.h"
//#include "bsnlib/internal/complex.h"
//#include "bsnlib/internal/operator.h"

namespace bsnlib {

namespace internal {

/** Transpose template operator definition */
template <typename S, typename A> class Operator<S, Operation<S, A, A, OperationsType::M_TRANSPOSE>>;

/** Transpose operator traits */
template <typename S, typename A> struct traits<Operator<S, Operation<S, A, A, OperationsType::M_TRANSPOSE>>>{
    enum{
        rows = traits<A>::cols,
        cols = traits<A>::rows
    };
};

/** Tranpose operation implementation */
template <typename S, typename A> struct OperationImplementation<S, A, A, OperationsType::M_TRANSPOSE>{
    static S eval(const A& a, const A&, const unsigned int i, const unsigned int j){
        return a(j, i);
    }
};

/** Transpose operator implementation */
template <typename Scalar, typename A> class Operator<Scalar, Operation<Scalar, A, A, OperationsType::M_TRANSPOSE>>{
public:

    typedef Operation<Scalar, A, A, OperationsType::M_TRANSPOSE> _BaseOperation;
    typedef Operator<Scalar, _BaseOperation> BaseElement;
    typedef Scalar S;

    /** constructor */
    Operator(const A& matrix) : Matrix(matrix){
        static_assert(traits<Scalar>::rows==1 && traits<Scalar>::cols==1, "Dimension error");
        static_assert(traits<A>::rows>1 || traits<A>::cols>1, "Dimension error");
    }

    /** entries acces operator */
    const Scalar operator()(const unsigned int i, const unsigned int j) const {
        assert(i>=0 && j>=0 && i<int(traits<A>::cols) && j<int(traits<A>::rows));
        return complex_conjugate<Scalar>::compute(Matrix(j,i));
        // return Matrix(j, i);
    }

    /** Addition operation */
    template <typename C> Operator<S, Operation<S, BaseElement, C, OperationsType::M_SUM_M>>
    operator + (const C& other) const {
        return Operator<S, Operation<S, BaseElement, C, OperationsType::M_SUM_M>>((*this), other);
    }

    /** Mult operation */
    template <typename C> Operator<S, Operation<S, BaseElement, C,
    traits<C>::rows==1?(traits<C>::cols==1?OperationsType::M_MUL_S:OperationsType::M_MUL_M):
        (OperationsType::M_MUL_M)>>
    operator * (const C& other) const {
        return Operator<S, Operation<S, BaseElement, C,
                                    traits<C>::rows==1?(traits<C>::cols==1?OperationsType::M_MUL_S:OperationsType::M_MUL_M):
                                        (OperationsType::M_MUL_M)>>((*this), other);
    }

    /** Substraction operator */
    template <typename C> Operator<S, Operation<S, BaseElement, C, OperationsType::M_SUB_M>>
    operator - (const C& other) const {
        return Operator<S, Operation<S, BaseElement, C, OperationsType::M_SUB_M>>((*this), other);
    }

protected:

    /** Matrix reference */
    const A& Matrix;
};

} // internal

} //bsnlib

#endif // TRANSPOSE_H
