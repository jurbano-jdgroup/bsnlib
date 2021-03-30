
#ifndef GIVENS_H
#define GIVENS_H

//#include "bsnlib/internal/traits.h"
//#include "bsnlib/internal/decompose.h"
//#include "bsnlib/solve/solve.h"
//#include "bsnlib/solve/backsubstitution.h"
//#include "bsnlib/solve/forwardsubstitution.h"
//#include "bsnlib/solve/directsubstitution.h"
//#include "bsnlib/internal/smatrix.h"
//#include "bsnlib/internal/slice.h"
//#include <math.h>

namespace bsnlib {

namespace solve {

/** Givens template definition */
template <typename S> struct Givens;

/** Givens based QR decomposition update */
template <typename S, typename _AType, typename _BType, int Rows, int Cols, int RIdx, int CIdx, bool hOder> struct QRGivensUpdate {};

/** Givens based QR decomposition update */
template <typename S, typename _AType, typename _BType, int Rows, int Cols, int RIdx, int CIdx> struct QRGivensUpdate<S,_AType,_BType,Rows,Cols,RIdx,CIdx,true>  {
    static bool compute(_AType& A, _BType& B) {
        S c, s;
        Givens<S>::givens(A(RIdx-2,CIdx-1), A(RIdx-1,CIdx-1), c, s);

        using internal::Slice;
        using internal::SMatrix;
        using internal::Operator;
        using internal::Operation;
        using internal::OperationsType;

        SMatrix<S,Rows,Rows> iden = SMatrix<S,Rows,Rows>::identity();
        iden(RIdx-2,RIdx-2)=c;
        iden(RIdx-2,RIdx-1)=s;
        iden(RIdx-1,RIdx-2)=-s;
        iden(RIdx-1,RIdx-1)=c;

        B *= iden;

        SMatrix<S,2,2> givensMatrix;
        givensMatrix(0,0)=c;
        givensMatrix(0,1)=s;
        givensMatrix(1,0)=-s;
        givensMatrix(1,1)=c;

        SMatrix<S,2,Cols-CIdx+1> tmp;
        Slice<S,_AType&,RIdx-2,RIdx-1,CIdx-1,Cols-1> slice(A);

        tmp = slice;

        using TOperator = Operator<S,Operation<S,SMatrix<S,2,2>,SMatrix<S,2,2>,OperationsType::M_TRANSPOSE>>;
        TOperator top(givensMatrix);

        using BOperator = Operator<S,Operation<S,TOperator,SMatrix<S,2,Cols-CIdx+1>,OperationsType::M_MUL_M>>;
        BOperator bop(top, tmp);

        slice = bop;

        constexpr int nRIdx = (RIdx-1) > CIdx ? (RIdx-1) : Rows;
        constexpr int nCIdx = (RIdx-1) > CIdx ? CIdx : (CIdx+1);
        constexpr bool cont = (Rows>Cols and nCIdx<=Cols) or (Rows==Cols && nCIdx<Cols);

        /* CIdx
         * x x x
         * x x x
         * x x x RIdx
         */

        //printf("RIdx:%d, CIdx:%d\n",RIdx,CIdx);

        return QRGivensUpdate<S,_AType,_BType,Rows,Cols,nRIdx,nCIdx,cont>::compute(A,B);
    }
};

/** Givens based QR decomposition update */
template <typename S, typename _AType, typename _BType, int Rows, int Cols, int RIdx, int CIdx> struct QRGivensUpdate<S,_AType,_BType,Rows,Cols,RIdx,CIdx,false>  {
    static bool compute(_AType&, _BType&) {
        return true;
    }
};

/** Givens rotations */
template <typename S> struct Givens {

    /** calculate the givens coefficients */
    static void givens(const S& a, const S& b, S& c, S& s) {
        if (b==static_cast<S>(0)) {
            c = static_cast<S>(1);
            s = static_cast<S>(0);
        }else{
            S tau = static_cast<S>(0);
            if (std::abs(b) > std::abs(a)) {
                tau = -a/b;
                s = 1.0 / std::sqrt(1.0 + tau*tau);
                c = s*tau;
            } else {
                tau = -b/a;
                c = 1.0 / std::sqrt(1.0 + tau*tau);
                s = c*tau;
            }
        }
    }

    /** calculate the givens matrix for the i, k indexes */
    template <typename _AMatrix, typename _BMatrix> static void matrix (const _AMatrix& A, _BMatrix& B, unsigned int i, unsigned int k, unsigned int j) {
        // assert(k==(i+1));
        using internal::traits;
        assert(i<int(traits<_AMatrix>::rows) && k<int(traits<_AMatrix>::rows) && j<int(traits<_AMatrix>::cols) && i<k);
        assert(i<int(traits<_BMatrix>::rows) && k<int(traits<_BMatrix>::rows) && i<int(traits<_BMatrix>::cols) && k<int(traits<_BMatrix>::cols) && i<k);

        S c, s;
        Givens<S>::givens(A(i,j), A(k,j), c, s);

        B(i, i) = c;
        B(i, k) = s;
        B(k, i) = -s;
        B(k, k) = c;
    }

    /** Givens rotations based QR decomposition */
    template <typename _MType,  int Rows, int Cols> static bool QRDecomposition(const _MType& A, internal::SMatrix<S, Rows, Rows>& Q, internal::SMatrix<S, Rows, Cols>& R) {
        /* check the dimensions */
        static_assert(Rows>2 && Cols>2 && Rows>=Cols, "Dimension error");
        /* check the matrix type */
        static_assert(int(internal::traits<_MType>::rows)==Rows && int(internal::traits<_MType>::cols)==Cols, "Dimension error");

        using internal::SMatrix;

        R = A;
        SMatrix<S,Rows,Rows>::makeIdentity(Q);

        return QRGivensUpdate<S,SMatrix<S,Rows,Cols>,SMatrix<S,Rows,Rows>,Rows,Cols,Rows,1,((Rows-2)>0)>::compute(R,Q);
    }
};

} // solve

} // bsnlib

#endif // GIVENS_H
