#ifndef RAW_H
#define RAW_H

//#include "bsnlib/internal/traits.h"
#include <assert.h>

namespace bsnlib {

namespace internal {

/** Raw operations */
struct RawOperation {
    enum {
        ROW_MAX                 = 1,
        COL_MAX                 = 2,
        COL_MAX_ABS             = 3,
        DIAGONAL_SUM            = 4,
        LOWER_SUM               = 5,
        OFF_DIAGONAL_SUM        = 6,
        OFF_DIAGONAL_SQR_SUM    = 7
    };
};

/** Raw operator definition */
template <typename S, int _Operation> struct Raw {};

/** Col max operator  */
template <typename S> struct Raw<S, RawOperation::COL_MAX> {

    template <typename _Matrix> static S compute(const _Matrix& A, unsigned int row, unsigned int col) {

        /* check the dimension at compilation time */
        static_assert (internal::traits<S>::rows==internal::traits<S>::cols, "Dimension error");
        /* check the matrix dimension at compilation time */
        static_assert (internal::traits<_Matrix>::rows>1 && internal::traits<_Matrix>::cols>1, "Dimension error");
        /* check the row and col specified */
        assert(col<int(internal::traits<_Matrix>::cols) && row<int(internal::traits<_Matrix>::rows));

        /* compute the col max */
        S tmp=static_cast<S>(0);
        for(unsigned int i=row; i<int(internal::traits<_Matrix>::rows); ++i){
            if(tmp < A(i, col)){
                tmp = A(i, col);
            }
        }

        return tmp;
    }
};

/** Col max abs operator  */
template <typename S> struct Raw<S, RawOperation::COL_MAX_ABS> {

    template <typename _Matrix> static S compute(const _Matrix& A, unsigned int row, unsigned int col) {

        /* check the dimension at compilation time */
        static_assert (internal::traits<S>::rows==internal::traits<S>::cols, "Dimension error");
        /* check the matrix dimension at compilation time */
        static_assert (internal::traits<_Matrix>::rows>1 && internal::traits<_Matrix>::cols>1, "Dimension error");
        /* check the row and col specified */
        assert(col<int(internal::traits<_Matrix>::cols) && row<int(internal::traits<_Matrix>::rows));

        /* compute the col max abs */
        S tmp=static_cast<S>(0);
        double abs_value=0.0;
        for(unsigned int i=row; i<int(internal::traits<_Matrix>::rows); ++i){
            abs_value=std::abs(A(i, col));
            if(tmp < abs_value){
                tmp = A(i, col);
            }
        }

        return tmp;
    }
};

/** Diagonal sum */
template <typename S> struct Raw<S, RawOperation::DIAGONAL_SUM> {

    template <typename _MatrixType> static S compute(const _MatrixType& A) {

        constexpr int rows = int(traits<_MatrixType>::rows);
        constexpr int cols = int(traits<_MatrixType>::cols);
        constexpr int dim = rows<cols?rows:cols;

        S tmp = static_cast<S>(0);

        for(unsigned int i=0; i<dim; ++i){
            tmp += A(i,i);
        }

        return tmp;
    }
};

/** Diagonal sum */
template <typename S> struct Raw<S, RawOperation::LOWER_SUM> {
    //     j
    //   x x x
    // i x x x
    //   x x x
    template <typename _MatrixType> static S compute(const _MatrixType& A) {
        constexpr int rows = int(traits<_MatrixType>::rows);
        constexpr int cols = int(traits<_MatrixType>::cols);

        S tmp = static_cast<S>(0);

        for(unsigned int i=0; i<rows; ++i){
            for(unsigned int j=0; j<cols; ++j){
                tmp += j<i?A(i,j):static_cast<S>(0);
            }
        }

        return tmp;
    }
};

/** Off diagonal sum */
template <typename S> struct Raw<S, RawOperation::OFF_DIAGONAL_SUM> {

    template <typename _MatrixType> static S compute (const _MatrixType &A) {

        /* check S */
        static_assert (is_scalar<S>::value, "S is not a scalar type");

        constexpr int rows = int(traits<_MatrixType>::rows);
        constexpr int cols = int(traits<_MatrixType>::cols);

        S tmp = static_cast<S>(0);

        for(unsigned int i=0; i<rows; ++i) {
            for (unsigned int j=0; j<cols; ++j) {
                tmp += i==j ? static_cast<S>(0) : A(i, j);
            }
        }

        return tmp;
    }
};

/** Off diagonal sum of pow 2 values */
template <typename S> struct Raw<S, RawOperation::OFF_DIAGONAL_SQR_SUM> {

    template <typename _MatrixType> static S compute (const _MatrixType& A) {
        // check S
        static_assert (is_scalar<S>::value, "S is not a scalar type");
        // check _MatrixType
        constexpr int rows = int(traits<_MatrixType>::rows);
        constexpr int cols = int(traits<_MatrixType>::cols);

        S tmp = static_cast<S>(0);
        using std::abs;

        for(unsigned int i=0; i<rows; ++i) {
            for(unsigned int j=0; j<cols; ++j) {
                tmp += i==j? static_cast<S>(0) : (abs(A(i,j)) * abs(A(i,j)));
            }
        }

        return tmp;
    }
};

} // internal

} // bsnlib

#endif // RAW_H
