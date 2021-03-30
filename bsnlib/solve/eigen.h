#ifndef EIGEN_H
#define EIGEN_H

/*#include "bsnlib/internal/traits.h"
#include "bsnlib/internal/smatrix.h"
#include "bsnlib/solve/householder.h"
#include "bsnlib/internal/raw.h"*/
#include <math.h>

namespace bsnlib {

namespace solve {

/** Eigen computation method */
struct EigenCompMethod {
    enum{
        DIMENSION_2=1,
        QR_ITERATION=2
    };
};

// max iteration in eigen computation
#define BSN_EIGEN_MAX_ITERS 100
// eigen computation tolerance
# define BSN_EIGEN_EPS 1e-33

/** Eigen computation solver */
template <typename S, int Method> struct EigenSolver {};

/** Second order dimension square matrix */
template <typename S> struct EigenSolver<S, EigenCompMethod::DIMENSION_2> {

    // Vector must be a complex values vector
    template <typename _Matrix, typename _Vector> static bool compute(const _Matrix& A, _Vector& v) {

        using internal::traits;
        constexpr int rows = traits<_Matrix>::rows;
        constexpr int cols = traits<_Matrix>::cols;
        constexpr int v_rows = traits<_Vector>::rows;

        static_assert (rows==cols && rows==2, "Matrix must be square");
        static_assert (rows==v_rows, "Vector rows must be the same as matrix rows");

        const S a = A(0,0);
        const S b = A(0,1);
        const S c = A(1,0);
        const S d = A(1,1);

        const S imag_sqr = (a+d)*(a+d) - 4.0*(a*d - b*c);
        const S real = (a+d)/2.0;

        S imag = static_cast<S>(0);

        if(imag_sqr < static_cast<S>(0)) {
            imag = std::sqrt((-1)*imag_sqr) / 2.0;

            v(0,0) = std::complex<S>(real, imag);
            v(1,0) = std::complex<S>(real, -1*imag);
        } else {
            imag = std::sqrt(imag_sqr) / 2.0;

            v(0,0) = real+imag;
            v(1,0) = real-imag;
        }

        return true;
    }
};

/** QR iteration eigen computation */
template <typename S> struct EigenSolver<S,EigenCompMethod::QR_ITERATION> {

    /** find p index */
    template <typename _Matrix> static void findPIndex(const _Matrix& H, int startIndex, int endIndex, int& p) {
        p = startIndex;

        while(p <= endIndex) {
            if(H(p+1,p) == static_cast<S>(0)) {
                p=p+1;
            } 
            else if(H(p+2,p+1) == static_cast<S>(0)) {
                p = p+2;
            } else {
                break;
            }
        }
    }

    /** find q index */
    template <typename _Matrix> static void findQIndex(const _Matrix& H, int startIndex, int endIndex, int& q) {
        q = endIndex;

        while(q > startIndex) {
            if(H(q, q-1) == static_cast<S>(0)) {
                q=q-1;
            }
            else if(H(q-1,q-2) == static_cast<S>(0)) {
                q=q-2;
            } else {
                break;
            }
        }
    }

    /** look for subdiagonal entries to set zero */
    template <typename _Matrix> static void cleanSubdiagonalEntries(_Matrix& H) {
        using internal::traits;
        constexpr int rows = traits<_Matrix>::rows;

        const S thres = static_cast<S>(1e-12);
        S eps = static_cast<S>(0);

        int i=0;
        while(i<(rows-2)) {
            eps = thres * ( std::abs(H(i,i)) + std::abs(H(i+1, i+1)) );
            if (std::abs(H(i+1,i)) < eps) {
                H(i+1, i) = static_cast<S>(0);
            } else if(std::abs(H(i+2,i+1)) < eps) {
                H(i+2,i+1) = static_cast<S>(0);
            }

            ++i;
        }
    }

    /** Compute the eigenvalues of a 2x2 matrix, Vector must have complex values entries */
    template <typename _Vector> static void compute2(_Vector& v, S a, S b, S c, S d) {
        const S imag_sqr = (a+d)*(a+d) - 4.0*(a*d - b*c);
        const S real = (a+d)/2.0;

        S imag = static_cast<S>(0);

        if(imag_sqr < static_cast<S>(0)) {
            imag = std::sqrt((-1)*imag_sqr) / 2.0;

            v(0,0) = std::complex<S>(real, imag);
            v(1,0) = std::complex<S>(real, -1*imag);
        } else {
            imag = std::sqrt(imag_sqr) / 2.0;

            v(0,0) = real+imag;
            v(1,0) = real-imag;
        }
    }

    /** compute the matrix eigen values, V must have complex entries */
    template <typename _Matrix, typename _Vector> static bool compute(const _Matrix& A, _Vector& V) {

        using _Solver = EigenSolver<S, EigenCompMethod::QR_ITERATION>;

        using internal::traits;
        constexpr int matrix_rows = int(traits<_Matrix>::rows);
        constexpr int matrix_cols = int(traits<_Matrix>::cols);

        // check square
        static_assert(matrix_cols==matrix_rows, "Matrix is not square");

        constexpr int v_rows = int(traits<_Vector>::rows);
        static_assert(matrix_rows<=v_rows, "Vector must have same rows as Matrix");
        static_assert(matrix_rows>2, "Greater than 2x2 matrix");

        // TODO: balance the matrix

        // compute the hessemberg reduction of the matrix
        using internal::SMatrix;
        using SMat = SMatrix<S, matrix_rows, matrix_rows>;
        SMatrix<S, matrix_rows, matrix_rows> H;

        if (!Householder<S>::template Hessemberg<_Matrix, SMatrix<S, matrix_rows, matrix_rows>>(A, H)) {
            return false;
        }

        SMat P0 = SMat::identity();

        using internal::Slice;
        Slice<S,SMat&,0,2,0,2> slice(P0);
        int n=matrix_rows-1, m=matrix_rows-2, s_i=0, k=0;

        int p=0;
        _Solver::template cleanSubdiagonalEntries<_Matrix>(H);
        _Solver::template findPIndex<_Matrix>(H, 0, matrix_rows-3, p);

        while(p<=(matrix_rows-3)) {

            S s, t, x, y, z;
            s = H(m,m) + H(n,n);
            t = H(m,m)*H(n,n) - H(m,n)*H(n,m);
            x = H(s_i,s_i)*H(s_i,s_i) + H(s_i,s_i+1)*H(s_i+1,s_i) - s*H(s_i,s_i) + t;
            y = H(s_i+1,s_i)*( H(s_i,s_i) + H(s_i+1,s_i+1) - s );
            z = H(s_i+1,s_i)*H(s_i+2,s_i+1);

            SMatrix<S,3,1> houseVector(x, y, z);
            SMatrix<S,3,3> houseMatrix = Householder<double>::template matrix<SMatrix<S,3,1>,3>(houseVector);

            slice = houseMatrix;
            SMat tmp = P0.transpose() * H * P0;

            //tmp.dump();
            
            _Solver::template cleanSubdiagonalEntries<SMat>(tmp);
            _Solver::template findPIndex<SMat>(tmp, 0, matrix_rows-3, p);
            //printf(", p index: %d\n", p);

            if(p>(matrix_rows-3)) {
                H = tmp;
                break;
            }

            if(!Householder<double>::template Hessemberg<SMat,SMat>(tmp, H)) {
                return false;
            }
            SMat::makeIdentity(P0);

            ++k;
        }

        //H.dump();

        // get the eigenvalues
        int i=0;
        while(i<(matrix_rows-2)) {
            if(H(i+1,i) == static_cast<S>(0)) {
                V(i,0) = H(i,i);
            } else if(H(i+2,i+1) == static_cast<S>(0)) {
                SMatrix<std::complex<S>,2,1> tmpVector;
                _Solver::template compute2<SMatrix<std::complex<S>,2,1>>(tmpVector, H(i,i), H(i,i+1),H(i+1,i),H(i+1,i+1));
                V(i,0) = tmpVector(0,0);
                V(i+1,0) = tmpVector(1,0);
                ++i;
            } else {
                return false;
            }

            ++i;
        }

        // check last indices
        if(i==n) {
            V(n,0) = H(n,n);
        } else if(i==m) {
            SMatrix<std::complex<S>,2,1> tmpVector;
            _Solver::template compute2<SMatrix<std::complex<S>,2,1>>(tmpVector, H(m,m), H(m,n), H(n,m), H(n,n));
            V(m,0) = tmpVector(0,0);
            V(n,0) = tmpVector(1,0);
        } else {
            return false;
        }

        return true;
    }
};

/** double QR iteration eigen solver */
typedef EigenSolver<double, EigenCompMethod::QR_ITERATION> EigenQrd;
typedef EigenSolver<double, EigenCompMethod::DIMENSION_2> Eigen2d;

} // solve

} // bsnlib

#endif // EIGEN_H