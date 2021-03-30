
#ifndef JACOBI_H
#define JACOBI_H

#include <assert.h>
#include <limits>
#include <cmath>
/*#include "bsnlib/internal/traits.h"
#include "bsnlib/internal/operator.h"
#include "bsnlib/internal/smatrix.h"*/

namespace bsnlib {

namespace solve {

/** Jacobi handler template struct */
template <typename S> struct Jacobi;

/** double jacobi specialization */
typedef Jacobi<double> Jacobid;

/** float jacobi specialization */
typedef Jacobi<float> Jacobif;

/** Jacobi implementation */
template <typename S> struct Jacobi {
    //===========================================
    /** symmetric rotation coefficients computation */
    template <typename _MatrixType> static void symmetric_rotation_coeffs(const _MatrixType &A, S &c, S &s, unsigned int p, unsigned int q) {
        // call using
        using internal::traits;
        using internal::is_scalar;
        /* check S */
        static_assert (is_scalar<S>::value, "S type is not a scalar type");
        /* check A */
        constexpr unsigned int a_rows = int(traits<_MatrixType>::rows);
        constexpr unsigned int a_cols = int(traits<_MatrixType>::cols);
        static_assert (a_rows>=2 && a_cols>=2 && a_rows>=a_cols, "_MatrixType is not a matrix type");
        /* check p && q */
        assert(p<a_cols && q<a_cols && p<q);

        const S u = A(p,q) - A(q,p);
        const S v = A(p,p) + A(q,q);

        // printf("u: %.3f, v: %.3f\n", u, v);

        if(std::numeric_limits<S>::min() > std::abs(u)) {
            s = static_cast<S>(0);
            c = static_cast<S>(1);
        } else {
            const S div = u/v;

            const S den = std::sqrt(static_cast<S>(1) + (div*div));
            c = static_cast<S>(1) / den;
            s = static_cast<S>(-1) * div * c;
        }
    }
    //=========================================
    /** rotation coefficients computation */
    template <typename _MatrixType> static void rotation_coeffs(const _MatrixType &A, S &c, S &s, unsigned int p, unsigned int q) {
        // call using
        using internal::traits;
        using internal::is_scalar;
        /* check S */
        static_assert (is_scalar<S>::value, "S type is not a scalar type");
        /* check A */
        constexpr unsigned int a_rows = int(traits<_MatrixType>::rows);
        constexpr unsigned int a_cols = int(traits<_MatrixType>::cols);
        static_assert (a_rows>=2 && a_cols>=2 && a_rows>=a_cols, "_MatrixType is not a matrix type");
        /* check p && q */
        assert(p<a_rows && q<a_cols && p<q);

        if(std::numeric_limits<S>::min() > std::abs(A(p,q))) {
            c = static_cast<S>(1);
            s = static_cast<S>(0);
        } else {
            const S tau = (A(q,q) - A(p,p)) / (2*A(p, q));
            S t;

            if (tau >= static_cast<S>(0)) {
                t = static_cast<S>(1) / (tau + std::sqrt(static_cast<S>(1) + tau*tau));
            }else {
                t = static_cast<S>(1) / (tau - std::sqrt(static_cast<S>(1) + tau*tau));
            }

            c = static_cast<S>(1) / std::sqrt(static_cast<S>(1) + t*t);
            s = t * c;
        }
    }
    //=========================================
    /** jacobi symmetric rotation matrix */
    template <typename _MatrixTypeA, typename _MatrixTypeB> static void symmetric_rotation_matrix (const _MatrixTypeA &A, _MatrixTypeB &B, unsigned int p, unsigned int q) {
        // call using
        using internal::is_scalar;
        using internal::traits;
        /* check S */
        static_assert (is_scalar<S>::value, "S is not a scalar type");
        /* check MA && MB */
        constexpr int a_rows = int(traits<_MatrixTypeA>::rows);
        constexpr int a_cols = int(traits<_MatrixTypeA>::cols);
        constexpr int b_rows = int(traits<_MatrixTypeB>::rows);
        constexpr int b_cols = int(traits<_MatrixTypeB>::cols);
        static_assert (a_rows>=2 && a_cols>=2 && a_rows>=a_cols && b_rows==2 && b_cols==2, "Matrix type A or B is not a matrix type");
        /* check q & p */
        assert(p<a_cols && q<a_cols && p<q);
        // Matrix type B must be an identity matrix
        S c, s;
        Jacobi<S>::template symmetric_rotation_coeffs<_MatrixTypeA>(A, c, s, p, q);
        // update B matrix
        B(0, 0) = c;
        B(0, 1) = s;
        B(1, 0) = -s;
        B(1, 1) = c;
    }
	//================================================
    /** jacobi rotation matrix computation */
    template <typename _MatrixTypeA, typename _MatrixTypeB> static void rotation_matrix (const _MatrixTypeA &A, _MatrixTypeB &B, unsigned int p, unsigned int q) {
        // call using
        using internal::is_scalar;
        using internal::traits;
        /* check S */
        static_assert (is_scalar<S>::value, "S is not a scalar type");
        /* check MA && MB */
        constexpr int a_rows = int(traits<_MatrixTypeA>::rows);
        constexpr int a_cols = int(traits<_MatrixTypeA>::cols);
        constexpr int b_rows = int(traits<_MatrixTypeB>::rows);
        constexpr int b_cols = int(traits<_MatrixTypeB>::cols);
        static_assert (a_rows>=2 && a_cols>=2 && a_rows>=a_cols && b_rows==2 && b_cols==2, "Matrix type A or B is not a matrix type");
        /* check q & p */
        assert(p<a_cols && q<a_cols && p<q);
        // Matrix type B must be an identity matrix
        S c, s;
        Jacobi<S>::template rotation_coeffs<_MatrixTypeA>(A, c, s, p, q);
        // update B matrix
        B(0, 0) = c;
        B(0, 1) = s;
        B(1, 0) = -s;
        B(1, 1) = c;
    }
    //================================================
    /** left and right rotation matrices */
    template <typename _MatrixTypeA, typename _MatrixTypeB> static void pair_matrices (const _MatrixTypeA &A, _MatrixTypeB &rot_left, _MatrixTypeB &rot_right, unsigned int p, unsigned int q) {
        // call using
        using internal::is_scalar;
        using internal::traits;
        using internal::SMatrix;
        /* check S */
        static_assert (is_scalar<S>::value, "S type is not a scalar type");
        /* check Matrix Type A */
        constexpr int a_rows = int(traits<_MatrixTypeA>::rows);
        constexpr int a_cols = int(traits<_MatrixTypeA>::cols);
        static_assert (a_rows>=2 && a_cols>=2 && a_rows>=a_cols, "Matrix Type A dimension error");
        /* check Matrix type B */
        constexpr int b_rows = int(traits<_MatrixTypeB>::rows);
        constexpr int b_cols = int(traits<_MatrixTypeB>::cols);
        static_assert (b_rows==2 && b_cols==2, "Matrix Type B is not a 2x2 matrix");
        /* check p & q */
        assert(p<a_cols && q<a_cols && p<q);
        // create aux matrices
        SMatrix<S,2,2> rot1, Am;
        Am(0,0)=A(p,p); Am(0,1)=A(p,q); Am(1,0)=A(q,p); Am(1,1)=A(q,q);
        Jacobi<S>::template symmetric_rotation_matrix<_MatrixTypeA, SMatrix<S,2,2>>(A, rot1, p, q);

        using MatrixBase = SMatrix<S,2,2>;
        using OpBase = internal::Operator<S, internal::Operation<S, MatrixBase, MatrixBase, internal::OperationsType::M_MUL_M>>;
        OpBase op(rot1, Am);

        Jacobi<S>::template rotation_matrix<OpBase, SMatrix<S,2,2>>(op, rot_right, 0, 1);

        rot_left = rot1.transpose() * rot_right;
    }
    //===============================================
    /** 2x2 jacobi svd solver */
    template <typename _MatrixTypeA, typename _MatrixTypeB> static void svd_2x2 (const _MatrixTypeA &A, _MatrixTypeB &E, _MatrixTypeB &U, _MatrixTypeB &V) {

        using internal::is_scalar;
        using internal::traits;
        using internal::SMatrix;

        /* check S */
        static_assert (is_scalar<S>::value, "S type is not a scalar type");
        /* check Matrix Type A */
        constexpr int a_rows = int(traits<_MatrixTypeA>::rows);
        constexpr int a_cols = int(traits<_MatrixTypeA>::cols);
        static_assert (a_rows==2 && a_cols==2, "Matrix Type A is not a 2x2 matrix");
        /* check Matrix type B */
        constexpr int b_rows = int(traits<_MatrixTypeB>::rows);
        constexpr int b_cols = int(traits<_MatrixTypeB>::cols);
        static_assert (b_rows==2 && b_cols==2, "Matrix Type B is not a 2x2 matrix");

        // A -> E
        // E = A;
        // I2x2 -> U
        U(0,0) = static_cast<S>(1); U(1,1) = static_cast<S>(1);
        U(0,1) = static_cast<S>(0); U(1,0) = static_cast<S>(0);
        // I2x2->V
        V(0,0) = static_cast<S>(1); V(1,1) = static_cast<S>(1);
        V(0,1) = static_cast<S>(0); V(1,0) = static_cast<S>(0);

        SMatrix<S,2,2> rot1, rot_right;
        Jacobi<S>::template symmetric_rotation_matrix<_MatrixTypeA, SMatrix<S,2,2>>(A, rot1, 0, 1);

        using MatrixBase = SMatrix<S,2,2>;
        using OpBase = internal::Operator<S, internal::Operation<S, MatrixBase, MatrixBase, internal::OperationsType::M_MUL_M>>;
        OpBase op(rot1, A);

        Jacobi<S>::template rotation_matrix<OpBase, SMatrix<S,2,2>>(op, rot_right, 0, 1);

        SMatrix<S,2,2> rot_left = rot1.transpose() * rot_right;

        E = rot_left.transpose() * A * rot_right;
        U *= rot_left;
        V *= rot_right;

        // Order of the singular values and the sign
        // are corrected in the general implementation
        // algorithm, so this method can not be used
        // to compute the 2x2 real svd matrix decom...
    }
};

} // solve

} // bsnlib

#endif // JACOBI_H
