
#ifndef HOUSEHOLDER_H
#define HOUSEHOLDER_H

/*#include "bsnlib/internal/traits.h"
#include "bsnlib/internal/dot.h"
#include "bsnlib/internal/slice.h"
#include "bsnlib/internal/smatrix.h"*/

namespace bsnlib{

namespace solve{

/** Householder struct definition */
template <typename S> struct Householder;

/** Double householder */
typedef Householder<double> Householderd;

/** Float householder */
typedef Householder<float> Householderf;

/** Householder QR decomposition update */
template <typename S, typename _AType, typename _VType, bool hOrder> struct UpdateQR {};

/** Householder QR decomposition update */
template <typename S, typename _AType, typename _VType> struct UpdateQR<S, _AType, _VType, true> {
	
    static bool compute(_AType& A, _VType& Q) {
		
        using internal::traits;
        constexpr unsigned int rows = static_cast<unsigned int>(traits<_AType>::rows);
        constexpr unsigned int cols = static_cast<unsigned int>(traits<_AType>::cols);
        constexpr unsigned int qRows = static_cast<unsigned int>(traits<_VType>::rows);

        using internal::Slice;
        using internal::SMatrix;

        Slice<S,const _AType&,0,rows-1,0,0> xVector(A);
        SMatrix<S,rows,1> vVector;
        S beta = static_cast<S>(0);

        Householder<S>::template vector<Slice<S,const _AType&,0,rows-1,0,0>, SMatrix<S,rows,1>>(xVector, vVector, beta);
        SMatrix<S,rows,rows> hMatrix(SMatrix<S,rows,rows>::identity() - (vVector*vVector.transpose()*beta));

        SMatrix<S,rows,cols> aMatrix;
        aMatrix = A;
        A = hMatrix * aMatrix;

        SMatrix<S,qRows,qRows> qUpdate = SMatrix<S,qRows,qRows>::identity();
        qUpdate.template slice<qRows-rows,qRows-1,qRows-rows,qRows-1>() = hMatrix;
        Q *= qUpdate;

        Slice<S,_AType&,1,rows-1,1,cols-1> aSlice(A);
        return UpdateQR<S, Slice<S,_AType&,1,rows-1,1,cols-1>, _VType, ((rows-1)>2)>::compute(aSlice, Q);
    }
};

/** Householder QR decomposition update, 2nd order specialization */
template <typename S, typename _AType, typename _VType> struct UpdateQR<S, _AType, _VType, false> {
	
    static bool compute(_AType& A, _VType& Q) {
		
        using internal::traits;
        constexpr unsigned int rows = static_cast<unsigned int>(traits<_AType>::rows);
        constexpr unsigned int cols = static_cast<unsigned int>(traits<_AType>::cols);
        constexpr unsigned int qRows = static_cast<unsigned int>(traits<_VType>::rows);

        using internal::Slice;
        using internal::SMatrix;

        Slice<S,const _AType&,0,rows-1,0,0> xVector(A);
        SMatrix<double,rows,1> vVector;
        S beta = static_cast<S>(0);

        Householder<S>::template vector<Slice<S,const _AType&,0,rows-1,0,0>, SMatrix<S,rows,1>>(xVector, vVector, beta);

        SMatrix<S,rows,rows> hMatrix(SMatrix<S,rows,rows>::identity() - (vVector*vVector.transpose()*beta));
        SMatrix<S,rows,cols> aMatrix;
        aMatrix = A;
        A = hMatrix * aMatrix;

        SMatrix<S,qRows,qRows> qUpdate = SMatrix<S,qRows,qRows>::identity();
        qUpdate.template slice<qRows-rows,qRows-1,qRows-rows,qRows-1>() = hMatrix;
        Q *= qUpdate;

        return true;
    }
};

/** Householder hessember reduction update */
template <typename S, typename _MatrixB, int index, bool hOrder> struct UpdateHessember {};

/** true hOrder specialization */
template <typename S, typename _MatrixB, int index> struct UpdateHessember<S, _MatrixB, index, true> {
    static bool compute(_MatrixB& H) {
        constexpr int rows = internal::traits<_MatrixB>::rows;
        using internal::SMatrix;
        using internal::Slice;
        constexpr int hMatrixRows = rows-index-1;

        using HouseSliceVector = Slice<S, _MatrixB&, index+1, rows-1, index, index>;
        HouseSliceVector slc(H);

        SMatrix<S, hMatrixRows, hMatrixRows> hMatrix = Householder<S>::template matrix<HouseSliceVector, hMatrixRows>(slc);

        Slice<S, _MatrixB&, index+1, rows-1, index, rows-1> sliceA(H);
        SMatrix<S, hMatrixRows, rows-index> matrixA;
        matrixA = sliceA;

        sliceA = hMatrix * matrixA;

        SMatrix<S, rows, hMatrixRows> hMat;
        Slice<S, _MatrixB&, 0, rows-1, index+1, rows-1> hMatSlc(H);
        hMat = hMatSlc;

        hMatSlc = hMat * hMatrix;

        constexpr int n_index = index+1;
        constexpr bool nHOrder = n_index<(rows-3);

        return UpdateHessember<S, _MatrixB, n_index, nHOrder>::compute(H);
    }
};

/** false hOrder specialization */
template <typename S, typename _MatrixB, int index> struct UpdateHessember<S, _MatrixB, index, false> {
    static bool compute(_MatrixB& H) {

        constexpr int rows = internal::traits<_MatrixB>::rows;
        using internal::SMatrix;
        using internal::Slice;
        constexpr int hMatrixRows = rows-index-1;

        using HouseSliceVector = Slice<S, _MatrixB&, index+1, rows-1, index, index>;
        HouseSliceVector slc(H);

        SMatrix<S, hMatrixRows, hMatrixRows> hMatrix = Householder<S>::template matrix<HouseSliceVector, hMatrixRows>(slc);

        Slice<S, _MatrixB&, index+1, rows-1, index, rows-1> sliceA(H);
        SMatrix<S, hMatrixRows, rows-index> matrixA;
        matrixA = sliceA;
        sliceA = hMatrix * matrixA;

        SMatrix<S, rows, hMatrixRows> hMat;
        Slice<S, _MatrixB&, 0, rows-1, index+1, rows-1> hMatSlc(H);
        hMat = hMatSlc;

        hMatSlc = hMat * hMatrix;

        return true;
    }
};

/** Householder matrix computation */
template <typename S> struct Householder {
	
    /** Get the householder vector v based on the x vector */
    template <typename _MType, typename _VType> static void vector(const _MType& x, _VType& v, S& betha) {

        using internal::traits;
        using internal::is_scalar;
        using internal::Dot;
        using internal::SMatrix;
        using std::sqrt;
        /* check the scalar type  */
        static_assert(is_scalar<S>::value, "Dimsension error");
        /* check the x vector  */
        static_assert (traits<_VType>::cols==1 && traits<_VType>::rows>1, "Dimension error");
        constexpr unsigned int R = static_cast<unsigned int>(traits<_VType>::rows);
        static_assert(int(internal::traits<_MType>::rows)==R && int(internal::traits<_MType>::cols)==1, "Dimension error");
        /* check the row dimension */        
        static_assert(R>1, "Dimension error");

        internal::Slice<S, const _MType&, 1, R-1, 0, 0> l(x);
        const S sigma = Dot<S>::template compute<internal::Slice<S, const _MType&, 1, R-1, 0, 0>, internal::Slice<S, const _MType&, 1, R-1, 0, 0>>(l, l);

        v(0,0) = static_cast<S>(1);

        for (unsigned int i=1; i<R; ++i) {
            v(i,0) = x(i,0);
        }

        if(sigma==static_cast<S>(0) && x(0, 0)>=static_cast<S>(0)) {
            betha = static_cast<S>(0);
        }
        else if(sigma==static_cast<S>(0) && x(0, 0)<static_cast<S>(0)) {
            betha = static_cast<S>(-2);
        }
        else {
            const S u = sqrt((x(0, 0)*x(0, 0)) + sigma);
            if(x(0, 0)<=static_cast<S>(0)) {
                v(0,0) = x(0, 0)- u;
            }
            else {
                v(0,0) = (-sigma) / (x(0, 0) + u);
            }

            betha = (static_cast<S>(2) * (v(0,0) * v(0,0))) / (sigma + (v(0,0) * v(0,0)));
            const S v_0_inv = static_cast<S>(1) / v(0,0);
            v *= v_0_inv;
        }
    }

    /** Get the householder matrix based on the x vector specified */
    template <typename _MType, int R> static internal::SMatrix<S, R, R> matrix(const _MType& x) {
        /* check the scalar type */
        static_assert(internal::traits<S>::rows==1 && internal::traits<S>::cols==1, "Dimension error");
        /* check the vector type */
        static_assert(int(internal::traits<_MType>::rows)==R && internal::traits<_MType>::cols==1, "Dimension error");
        /* check the vector dimension */
        static_assert(R>1, "Dimension error");
        internal::SMatrix<S,R,1> v;
        S beta;

        Householder::template vector<_MType, internal::SMatrix<S,R,1>>(x, v, beta);
        return internal::SMatrix<S,R,R>(internal::SMatrix<S,R,R>::identity() - (v*v.transpose()*beta));
    }


    /** Householder based QR decomposition */
    template <typename _MType, int Rows, int Cols> static bool QRDecomposition(const _MType& A, internal::SMatrix<S, Rows, Rows>& Q, internal::SMatrix<S, Rows, Cols>& R) {
        /* check the dimensions */
        static_assert(Rows>2 && Cols>2 && Rows>=Cols, "Dimension error");
        /* check the matrix type */
        static_assert(int(internal::traits<_MType>::rows)==Rows && int(internal::traits<_MType>::cols)==Cols, "Dimension error");

        internal::SMatrix<S, Rows, Rows>::makeIdentity(Q);
        R = A;

        return UpdateQR<S, internal::SMatrix<S,Rows,Cols>, internal::SMatrix<S,Rows,Rows>, (Rows>2)>::compute(R, Q);
    }

    /** Householder based hessember reduction */
    template <typename _MatrixA, typename _MatrixH> static bool Hessemberg(const _MatrixA& A, _MatrixH& H) {
        // check dimensions
        constexpr int a_rows = internal::traits<_MatrixA>::rows;
        constexpr int a_cols = internal::traits<_MatrixA>::cols;
        constexpr int h_rows = internal::traits<_MatrixH>::rows;
        constexpr int h_cols = internal::traits<_MatrixH>::cols;

        static_assert (a_rows==a_cols && a_rows>2, "A matrix must be square");
        static_assert (a_rows==h_rows && a_cols==h_cols, "A matrix and H matrix mus have the same dimension");

        H = A;
        constexpr bool hOrder = a_rows>3;

        return UpdateHessember<S, _MatrixH, 0, hOrder>::compute(H);
    }
};

} // solve

} // bsnlib

#endif // HOUSEHOLDER_H
