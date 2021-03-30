
#ifndef OPERATOR_H
#define OPERATOR_H

namespace bsnlib{

namespace internal{

/** operator */
template <typename Scalar, typename Operation> class Operator{};

/** operations type */
struct OperationsType {
    enum {
        M_SUM_M = 0,
        M_SUB_M = 1,
        M_MUL_M = 2,
        M_MUL_S = 3,
        M_TRANSPOSE = 4
    };
};

/** operation struct */
template <typename Scalar, typename A, typename B, int Op> struct Operation {};

/** operation implementation */
template <typename Scalar, typename A, typename B, int Op> struct OperationImplementation{};

} // internal

} // bsnlib

#endif // OPERATOR_H
