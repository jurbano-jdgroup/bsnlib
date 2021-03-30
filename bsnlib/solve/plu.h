
#ifndef PLU_H
#define PLU_H

//#include "bsnlib/internal/smatrix.h"
//#include "bsnlib/solve/solve.h"
//#include "bsnlib/internal/basern.h"
//#include "bsnlib/internal/slice.h"

namespace bsnlib {

namespace solve {

/** pivoting plu implementation type */
struct _PLU_Method_Type {
    enum {
        PARTIAL_PIVOTING=1,
        COMPLETE_PIVOTING=2
    };
};

/** PLU decomposition implementation */
template <typename S, int _Method> struct _PLU_Implementation {};

/** */
template <typename S> struct _PLU_Implementation<S, _PLU_Method_Type::PARTIAL_PIVOTING> {
    template <typename _Matrix, int R, int C> static bool compute(const _Matrix& A, internal::SMatrix<S, R, 1>& P, internal::SMatrix<S, R, C>&L, internal::SMatrix<S, R, C>& U ) {
        /* check S */
        static_assert(internal::traits<S>::rows==internal::traits<S>::cols, "Dimension error");
        /* check the dimensions of the matrices */
        static_assert(internal::traits<_Matrix>::rows==R && internal::traits<_Matrix>::cols==C, "Dimension error");
        /* check the allowed dimension */
        static_assert(R==C && C>=3, "Dimension error");

        /* compute the decomposition */
        U = A;
        L.setZero();
        P.setZero();

        unsigned int u=0;
        S row_max = static_cast<S>(0);
        S tmp_val = static_cast<S>(0);
        for(unsigned int k=0; k<R-1; ++k){
            /* find the index */
			row_max=U(k,k);
            u=k;
            for(unsigned int i=k+1; i<R; ++i){
                if(U(i, k)>row_max){
                    row_max = U(i, k);
                    u=i;
                }
            }
            /* update the pivot value */
            P(k,0)=u;
            /* switch the rows */
            U.switchRows(k,u);
            /* update the matrix */
            tmp_val=U(k,k);
            assert(tmp_val!=0);
            /* update the k left vector */
            for(unsigned int i=k+1; i<R; ++i){
                U(i,k) *= 1.0 / tmp_val;
            }
            /* update the next matrix */
            for(unsigned int i=k+1; i<R; ++i){
                for(unsigned int j=k+1; j<R; ++j){
                    U(i, j) -= U(i, k)*U(k, j);
                }
            }
        }

        /** update the matrices */
        for(unsigned int i=0;i<R;++i){
            for(unsigned int j=0; j<C; ++j){
                if(i==j){
                    L(i,j)=1.;
                }else if(i>j){
                    L(i,j)=U(i,j);
                    U(i,j)=0.;
                }
            }
        }

        return true;
    }
};

/** PLU inverse update */
template <typename S, int Rows, int Idx, bool IdxLeft> struct PluInvUpdate {};

/** PLU inverse update
 * Solve the system Ax=b
 * L(Ux) = Pb, Ux = y
 * Ly = Pb,
 * Ux = y
 */
template <typename S, int Rows, int Idx> struct PluInvUpdate<S, Rows, Idx, true> {
    template <typename _MatrixType, typename _VectorType> static void compute(_MatrixType& A, const _MatrixType& P, const _MatrixType& L, const _MatrixType& U, _VectorType& b){
        using _RealBase = internal::RnBase<S,Rows,Idx>;
        _RealBase e;
        using _BaseOperator = internal::Operator<S, internal::Operation<S, _MatrixType, _RealBase, internal::OperationsType::M_MUL_M>>;
        _BaseOperator op(P, e);

        internal::SMatrix<S,Rows,1> y;
        Solve<S, SolveMethod::FORWARD>::template solve<_MatrixType, internal::SMatrix<S,Rows,1>, _BaseOperator>(L, y, op);
        Solve<S, SolveMethod::BACK>::template solve<_MatrixType, _VectorType, internal::SMatrix<S,Rows,1>>(U, b, y);

        internal::Slice<S,_MatrixType&,0,Rows-1,Idx,Idx> slice(A);

        return PluInvUpdate<S,Rows,Idx+1,((Rows-Idx-1)>0)>::template compute<_MatrixType, internal::Slice<S,_MatrixType&,0,Rows-1,Idx,Idx>>(A, P, L, U, slice);
    }
};

/** PLU inverse update */
template <typename S, int Rows, int Idx> struct PluInvUpdate<S, Rows, Idx, false> {
    template <typename _MatrixType, typename _VectorType> static void compute(_MatrixType& A, const _MatrixType& P, const _MatrixType& L, const _MatrixType& U, _VectorType& b){
        using _RealBase = internal::RnBase<S,Rows,Idx>;
        _RealBase e;
        using _BaseOperator = internal::Operator<S, internal::Operation<S, _MatrixType, _RealBase, internal::OperationsType::M_MUL_M>>;
        _BaseOperator op(P, e);

        internal::SMatrix<S,Rows,1> y;
        Solve<S, SolveMethod::FORWARD>::template solve<_MatrixType, internal::SMatrix<S,Rows,1>, _BaseOperator>(L, y, op);
        Solve<S, SolveMethod::BACK>::template solve<_MatrixType, _VectorType, internal::SMatrix<S,Rows,1>>(U, b, y);

        return;
    }
};

/** double partial pivoting plu implementation */
typedef _PLU_Implementation<double, _PLU_Method_Type::PARTIAL_PIVOTING> PLUImpd;

/** double 3x3 matrix decomposition method */
typedef Decompose<double, DecompositionMethod::PLU, 3, 3> PLU3d;

/** double 3x4 matrix decomposition method */
typedef Decompose<double, DecompositionMethod::PLU, 4, 4> PLU4d;

/** float 3x3 matrix decomposition method */
typedef Decompose<float, DecompositionMethod::PLU, 3, 3> PLU3f;

/** float 3x4 matrix decomposition method */
typedef Decompose<float, DecompositionMethod::PLU, 4, 4> PLU4f;

/** PLU decomposition method */
template <typename S, int R, int C> class Decompose<S, DecompositionMethod::PLU, R, C> {
public:

    typedef internal::SMatrix<S,R,C> _BaseMatrix;
    typedef internal::SMatrix<S,R,1> _BaseVector;

    /** constructor */
    Decompose(const internal::SMatrix<S, R, C>& _A) : decomposed(false), A(_A), P("P"), L("L"), U("U") {
        /* check the scalar type */
        static_assert(internal::traits<S>::cols==1 && internal::traits<S>::rows==1, "Dimension error");
        /* check the dimensions */
        static_assert(R==C && C>=3, "Dimension error");
        /* do the decomposition */
        decomposed=doDecomposition();
    }

    /** do the matrix decomposition */
    inline bool doDecomposition() {
        return _PLU_Implementation<S, _PLU_Method_Type::PARTIAL_PIVOTING>::template compute<internal::SMatrix<S,R,C>,R,C>(A,P,L,U);
    }

    /** Get the L matrix */
    inline const internal::SMatrix<S,R,C>& getL() const {
        assert(decomposed);
        return L;
    }

    /** Get the U matrix */
    inline const internal::SMatrix<S,R,C>& getU() const {
        assert(decomposed);
        return U;
    }

    /** Get the P vector */
    inline const internal::SMatrix<S,R,1>& getPVector() const {
        assert(decomposed);
        return P;
    }

    /** Get the P matrix */
    inline internal::SMatrix<S,R,C> getP() const {
        assert(decomposed);
        internal::SMatrix<S,R,C> ret = internal::SMatrix<S,R,C>::identity();
        unsigned int r0=0, r1=0;
        for(unsigned int k=0; k<R-1; k++){
            r0 = k;
            r1 = P[k];
            ret.switchRows(r0<r1?r0:r1, r0<r1?r1:r0);
        }

        return ret;
    }

    /** Solve the system Ax=b
     * L(Ux) = Pb, Ux = y
     * Ly = Pb,
     * Ux = y
     */
    template <typename _VectorType> inline internal::SMatrix<S,R,1> solve(const _VectorType& b) const {
        /* check if decomposed */
        assert(decomposed && R==C && R>=3);
        using _BaseOperator = internal::Operator<S, internal::Operation<S, _BaseMatrix, _VectorType, internal::OperationsType::M_MUL_M>>;
        _BaseOperator op((*this).getP(), b);

        internal::SMatrix<S,R,1> y;
        Solve<S, SolveMethod::FORWARD>::template solve<_BaseMatrix, _BaseVector, _BaseOperator>(L, y, op);
        internal::SMatrix<S,R,1> x;
        Solve<S, SolveMethod::BACK>::template solve<_BaseMatrix, _BaseVector, _BaseVector>(U, x, y);
        return x;
    }

    /** Solve the system Ax=b */
    template <typename _VectorType, typename _VectorType2> void solve(const _VectorType& b, _VectorType2& x) const {
        /* check if decomposed */
        assert(decomposed && R==C && R>=3);
        using _BaseOperator = internal::Operator<S, internal::Operation<S, _BaseMatrix, _VectorType, internal::OperationsType::M_MUL_M>>;
        _BaseOperator op((*this).getP(), b);

        internal::SMatrix<S,R,1> y;
        Solve<S, SolveMethod::FORWARD>::template solve<_BaseMatrix, _BaseVector, _BaseOperator>(L, y, op);
        //internal::SMatrix<S,R,1> x;
        Solve<S, SolveMethod::BACK>::template solve<_BaseMatrix, _BaseVector, _BaseVector>(U, x, y);
    }

    /** inverse */
    inline internal::SMatrix<S, R, C> inverse() const {
        /* check if decomposed */
        assert(decomposed && R==C && R>=3);
        internal::SMatrix<S,R,C> ret;
        internal::Slice<S,internal::SMatrix<S,R,C>&,0,R-1,0,0> slice(ret);

        PluInvUpdate<S,R,1,((R-1)>0)>::template compute<internal::SMatrix<S,R,C>, internal::Slice<S,internal::SMatrix<S,R,C>&,0,R-1,0,0>>(ret,(*this).getP(),L,U,slice);

        return ret;
    }

    /** determinant */
    inline S determinant() const {
        static_assert (R==C, "Determinant on non square matrix");
        assert(decomposed);
        S tmp = static_cast<S>(1);
        for(unsigned k=0; k<R; ++k){
            tmp *= U(k,k);
        }
        for(unsigned k=0; k<R-1; ++k){
            tmp *= P(k,0)==static_cast<S>(k)?static_cast<S>(1):static_cast<S>(-1);
        }
        return tmp;
    }

    /** is decomposed */
    inline const bool& isDecomposed() const {
        return decomposed;
    }

protected:

    /** decomposed */
    bool decomposed;
    /** A matrix */
    const internal::SMatrix<S, R, C>& A;
    /** P matix */
    internal::SMatrix<S, R, 1> P;
    /** L matrix */
    internal::SMatrix<S, R, C> L;
    /** U matrix */
    internal::SMatrix<S, R, C> U;
};

} // solve

} // bsnlib

#endif // PLU_H
