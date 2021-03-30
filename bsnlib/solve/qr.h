#ifndef QR_H
#define QR_H

//#include "bsnlib/internal/decompose.h"
//#include "bsnlib/internal/operator.h"
//#include "bsnlib/internal/smatrix.h"
//#include "bsnlib/solve/solve.h"
//#include "bsnlib/internal/basern.h"

namespace bsnlib {

namespace solve {

/** QR decomposition method */
struct QR_METHOD {
    enum {
        HOUSEHOLDER = 0,
        GIVENS = 1
    };
};

/** QR inverse update */
template <typename S, int Rows, int Idx, bool IdxLeft> struct QRInvUpdate {};

/** QR inverse update */
template <typename S, int Rows, int Idx> struct QRInvUpdate<S, Rows, Idx, true> {
    /**
     * Solve the system Ax = b
     * QR x = b,
     * Rx = Qt b
     */
    template <typename _MatrixType, typename _VectorType> static void compute(_MatrixType& A, const _MatrixType& Q, const _MatrixType& R, _VectorType& b) {
        using _RnBase = internal::RnBase<S,Rows,Idx>;
        _RnBase e;
        using _TOperator = internal::Operator<S,internal::Operation<S,_MatrixType,_MatrixType,internal::OperationsType::M_TRANSPOSE>>;
        _TOperator top(Q);
        using _MOperator = internal::Operator<S,internal::Operation<S,_TOperator,_RnBase,internal::OperationsType::M_MUL_M>>;
        _MOperator mop(top,e);
        Solve<S,SolveMethod::BACK>::template solve<_MatrixType,_VectorType,_MOperator>(R,b,mop);

        internal::Slice<S,_MatrixType&,0,Rows-1,Idx,Idx> slice(A);

        return QRInvUpdate<S,Rows,Idx+1,((Rows-Idx-1)>0)>::template compute<_MatrixType, internal::Slice<S,_MatrixType&,0,Rows-1,Idx,Idx>>(A, Q, R, slice);
    }
};

/** QR inverse update */
template <typename S, int Rows, int Idx> struct QRInvUpdate<S, Rows, Idx, false> {
    /**
     * Solve the system Ax = b
     * QR x = b,
     * Rx = Qt b
     */
    template <typename _MatrixType, typename _VectorType> static void compute(_MatrixType& A, const _MatrixType& Q, const _MatrixType& R, _VectorType& b) {
        using _RnBase = internal::RnBase<S,Rows,Idx>;
        _RnBase e;
        using _TOperator = internal::Operator<S,internal::Operation<S,_MatrixType,_MatrixType,internal::OperationsType::M_TRANSPOSE>>;
        _TOperator top(Q);
        using _MOperator = internal::Operator<S,internal::Operation<S,_TOperator,_RnBase,internal::OperationsType::M_MUL_M>>;
        _MOperator mop(top,e);
        Solve<S,SolveMethod::BACK>::template solve<_MatrixType,_VectorType,_MOperator>(R,b,mop);

        return;
    }
};

/** QR decomposition definition */
template <typename S, int R, int C> class Decompose<S, DecompositionMethod::QR, R, C>;

/** double 3x3 matrix qr decomposition */
typedef Decompose<double, DecompositionMethod::QR, 3, 3> QR3d;

/** double 4x4 matrix qr decomposition */
typedef Decompose<double, DecompositionMethod::QR, 4, 4> QR4d;

/** float 3x3 matrix qr decomposition */
typedef Decompose<float, DecompositionMethod::QR, 3, 3> QR3f;

/** float 4x4 matrix qr decomposition */
typedef Decompose<float, DecompositionMethod::QR, 4, 4> QR4f;

/** QR matrix decomposition implementation */
template <typename S, int R, int C> class Decompose<S, DecompositionMethod::QR, R, C> {
public:

    typedef internal::SMatrix<S,R,C> _MatrixType;
    typedef internal::SMatrix<S,R,1> _VectorType;

    /** Default constructor */
    Decompose(const internal::SMatrix<S, R, C>& matrix, bool decompose=false) : M(matrix) {
        /* check the s type */
        static_assert(internal::traits<S>::cols==1 && internal::traits<S>::rows==1, "Dimension error");
        /* min 3x3 matrix dimensions allowed */
        static_assert(R>=C && C>2, "Dimension error");
        /* check if the decomposition must be perfomed */
        if(decompose){
            decomposed = doDecomposition<QR_METHOD::HOUSEHOLDER>();
        }
    }

    /** Solve the system Ax=y only if the matrix A is square */
    template <typename _OtherVectorType> inline internal::SMatrix<S,R,1> solve(const _OtherVectorType& y) {
        /* check the dimesion */
        assert(R==C && R>=3 && decomposed);

        using _TOperator = internal::Operator<S, internal::Operation<S, _MatrixType, _MatrixType, internal::OperationsType::M_TRANSPOSE>>;
        _TOperator top(QMatrix);
        using _BaseOperator = internal::Operator<S, internal::Operation<S, _TOperator, _OtherVectorType, internal::OperationsType::M_MUL_M>>;
        const _BaseOperator op(top, y);

        internal::SMatrix<S,R,1> x;
        Solve<S, SolveMethod::BACK>::template solve<_MatrixType, _VectorType, _BaseOperator>(RMatrix, x, op);
        return x;
    }

    /** Do the matrix decomposition and return true if success */
    template <int Method> inline bool doDecomposition() {
        if(Method==QR_METHOD::HOUSEHOLDER){
            return Householder<S>::template QRDecomposition<internal::SMatrix<S,R,C>, R,C>(M, QMatrix, RMatrix);
        }else if(Method==QR_METHOD::GIVENS){
            return Givens<S>::template QRDecomposition<internal::SMatrix<S,R,C>, R,C>(M, QMatrix, RMatrix);
        }else{
            return false;
        }
    }

    /** Get the Q Matrix */
    inline const internal::SMatrix<S,R,R>& getQ() const {
        assert(decomposed);
        return QMatrix;
    }

    /** Get the R Matrix */
    inline const internal::SMatrix<S,R,C> getR() const {
        assert(decomposed);
        return RMatrix;
    }

    /** inverse */
    inline internal::SMatrix<S,R,C> inverse() const {
        static_assert (R==C, "Inverse on non square matrix");
        assert(decomposed);

        internal::SMatrix<S,R,C> ret;
        internal::Slice<S,internal::SMatrix<S,R,C>&,0,R-1,0,0> slice(ret);

        QRInvUpdate<S,R,1,((R-1)>0)>::template compute<internal::SMatrix<S,R,C>, internal::Slice<S,internal::SMatrix<S,R,C>&,0,R-1,0,0>>(ret, QMatrix, RMatrix, slice);

        return ret;
    };

    /** determinant */
    inline S determinant () const {
        static_assert (R==C, "Determinant on non square matrix");
        assert(decomposed);

        S tmp = static_cast<S>(1);
        for(unsigned int k=0; k<R; ++k){
            tmp *= RMatrix(k,k);
        }

        return tmp;
    }

    /** Check if decomposed */
    inline bool isDecomposed() {
        return decomposed;
    }

protected:

    /** Matrix */
    const internal::SMatrix<S,R,C>& M;
    /** Orthogonal matrix */
    internal::SMatrix<S,R,R> QMatrix;
    /** Upper triangular matrix */
    internal::SMatrix<S,R,C> RMatrix;
    /** Decompose flag */
    bool decomposed = false;
};


} // solve

} // bsnlib

#endif // QR_H
