
#ifndef SVD_H
#define SVD_H

#include <limits>
/*#include "bsnlib/internal/traits.h"
#include "bsnlib/internal/raw.h"
#include "bsnlib/internal/smatrix.h"
#include "bsnlib/solve/jacobi.h"*/

namespace bsnlib {

namespace solve {

/** SVD Solvers */
struct SVD_Solver_Name {
    enum {
        Jacobi = 1
    };
};

/** SVD decomposition struct */
template <typename S, int Solver> struct SVD_Solver {};

/** double SVD solver specialization */
typedef SVD_Solver<double, SVD_Solver_Name::Jacobi> SVDSolverd;

/** float SVD solver specialization */
typedef SVD_Solver<float, SVD_Solver_Name::Jacobi> SVDSolverf;

/** Jacobi SVD Solver */
template <typename S> struct SVD_Solver<S, SVD_Solver_Name::Jacobi> {

    /** Jacobi based SVD decomposition computation */
    template <typename _MatrixTypeA, typename _MatrixTypeU, typename _MatrixTypeV> static void compute (_MatrixTypeA &A, _MatrixTypeU &U, _MatrixTypeV &V) {
        // check S
        using internal::is_scalar;
        static_assert (is_scalar<S>::value, "S is not a scalar type");
        // check _Matrix Type
        using internal::traits;
        constexpr int rows = int(traits<_MatrixTypeA>::rows);
        constexpr int cols = int(traits<_MatrixTypeA>::cols);
        static_assert (rows>=2 && cols>=2, "Matrix Type is not a matrix type");
        // check Matrix Type for U
        static_assert (int(traits<_MatrixTypeU>::rows)==rows && int(traits<_MatrixTypeU>::cols)==cols, "U dimensions error");
        // check Matrix Type for V
        static_assert (int(traits<_MatrixTypeV>::rows)==cols && int(traits<_MatrixTypeV>::cols)==cols, "V dimensions error");
        // square matrix for now
        static_assert (rows==cols, "Square matrices for now");

        // TODO: implement QR preconditioner

        using internal::Raw;
        using internal::RawOperation;
        S prec = static_cast<S>(std::numeric_limits<S>::min());
        S maxVal = Raw<S, RawOperation::OFF_DIAGONAL_SQR_SUM>::template compute<_MatrixTypeA>(A);
        S thres = prec * maxVal;
        // U && V must be identity
        using internal::SMatrix;
        using Mat2 = internal::SMatrix<S,2,2>;

        Mat2 rot_left, rot_right;
        SMatrix<S,rows,cols> Rleft, Rright, tmpB;

        SMatrix<S,rows,cols>::makeIdentity(Rleft);
        SMatrix<S,rows,cols>::makeIdentity(Rright);

        S entries_sum;
        bool finished = false;

        /**** iteration ****/
        while(!finished) {

            entries_sum = static_cast<S>(0);

            for(unsigned int p=0; p<(cols-1); ++p) {
                for(unsigned int q=(p+1); q<cols; ++q) {
                    entries_sum += (A(p,q)*A(p,q)) + (A(q,p)*A(q,p));
                    Jacobi<S>::template pair_matrices<_MatrixTypeA, Mat2>(A, rot_left, rot_right, p, q);
                    // update Rleft
                    Rleft(p,p)=rot_left(0,0);
                    Rleft(p,q)=rot_left(0,1);
                    Rleft(q,p)=rot_left(1,0);
                    Rleft(q,q)=rot_left(1,1);
                    // update Rright
                    Rright(p,p)=rot_right(0,0);
                    Rright(p,q)=rot_right(0,1);
                    Rright(q,p)=rot_right(1,0);
                    Rright(q,q)=rot_right(1,1);
                    // update B
                    tmpB = Rleft.transpose() * A * Rright;
                    A = tmpB;
                    // update U & V
                    U *= Rleft;
                    V *= Rright;
                    // clean Rleft & Rright
                    SMatrix<S,rows,cols>::makeIdentity(Rleft);
                    SMatrix<S,rows,cols>::makeIdentity(Rright);
                }
            }

            finished = entries_sum <= thres;
            //printf("sum value: %.7f, thres: %.7f\n", entries_sum, thres);
        }
        /*** sign correction ***/
        for(unsigned int i=0; i<cols; ++i){
            if(A(i,i) < static_cast<S>(0)){
                A(i,i) *= static_cast<S>(-1);
                for(unsigned int j=0; j<rows; ++j){
                    U(j,i) *= static_cast<S>(-1);
                }
            }
        }
        /*** order correction ***/
        unsigned int maxPos = 0;
        S tmpVal = static_cast<S>(0);
        for(unsigned int i=0; i<cols; ++i)
        {
            // search max pos
            maxPos = i;
            for(unsigned int j=(i+1); j<cols; ++j){
                if (A(j,j) > A(maxPos,maxPos)) {
                    maxPos = j;
                }
            }
            // update
            if (maxPos != i) {
                tmpVal = A(i,i);
                A(i,i) = A(maxPos, maxPos);
                A(maxPos, maxPos) = tmpVal;
                // update cols for U
                for (unsigned int j=0; j<rows; ++j){
                    tmpVal = U(j,i);
                    U(j,i)=U(j,maxPos);
                    U(j,maxPos)=tmpVal;
                }
                // update cols for V
                for(unsigned int j=0; j<cols; ++j){
                    tmpVal = V(j,i);
                    V(j,i)=V(j,maxPos);
                    V(j,maxPos)=tmpVal;
                }
            }
        }
    }
};

} // solve

} // bsnlib

#endif // SVD_H
