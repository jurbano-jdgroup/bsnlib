#ifndef LDLT_H
#define LDLT_H

//#include "bsnlib/internal/traits.h"
//#include "bsnlib/internal/decompose.h"
//#include "bsnlib/solve/solve.h"
//#include "bsnlib/internal/smatrix.h"

namespace bsnlib {

namespace solve {

/** LDLT method type */
struct _LDLT_Method_Type{
    enum{
        PARTIAL_PIVOTING = 1
    };
};

/** LDLT decomposition implementation */
template <typename S, int _Method> struct _LDLT_Implementation {};

/** */
template <typename S> struct _LDLT_Implementation<S, _LDLT_Method_Type::PARTIAL_PIVOTING>{
    /** computation */
    template <typename _MatrixType, int R, int C> static bool compute(const _MatrixType& A, internal::SMatrix<S,R,1>& P, internal::SMatrix<S,R,C>& L, internal::SMatrix<S,R,C>& D) {
        /* check scalar type */
        static_assert (internal::traits<S>::rows==1 && internal::traits<S>::cols==1, "Scalar dimension error");
        /* check matrix type */
        static_assert (internal::traits<_MatrixType>::rows==R && internal::traits<_MatrixType>::cols==C, "_MatrixType dimension error");
        static_assert (R==C && R>=3, "Must be a square matrix");

        L = A;
        P.setZero();
        D.setZero();

        S max_aux = static_cast<S>(0);
        S alpha = static_cast<S>(0);
        unsigned int index = 0;

        for(unsigned int k=0; k<R-1; ++k){
            /* find the max index */
            max_aux=L(k,k);
            index = 0;
            for(unsigned int i=k+1; i<R; ++i){
                if(L(i, i)>max_aux){
                    max_aux = L(i,i);
                    index = i;
                }
            }
            /* update the pivot */
            P(k,0)=index;
            /* switch row and col */
            L.switchRows(k,index);
            L.switchCols(k,index);
            /* get and check alpha */
            alpha = L(k,k);
            assert(alpha!=0);
            /* update row vector */
            for(unsigned int i=k+1; i<R; ++i){
                L(i,k) *= 1.0/alpha;
            }
            /* update the next matrix */
            for(unsigned int i=k+1; i<R; ++i){
                for(unsigned int j=k+1; j<R; ++j){
                    L(i, j) -= L(i, k)*L(k, j);
                }
            }
        }

        /* update the matrices */
        for(unsigned int i=0;i<R;++i){
            for(unsigned int j=0; j<C; ++j){
                if(i==j){
                    D(i,j)=L(i,j);
                    L(i,j)=1.;
                }else if(j>i){
                    L(i,j)=0.0;
                }
            }
        }

        return true;
    }
};

/** double partial pivoting ldlt implementation */
typedef _LDLT_Implementation<double, _LDLT_Method_Type::PARTIAL_PIVOTING> LDLTImpd;

/** LDLT decomposition class definition */
template <typename S, int R, int C> class Decompose<S, DecompositionMethod::LDLT, R, C>;

/** LDLT, 3x3 float matrix decomposition short definition */
typedef Decompose<float, DecompositionMethod::LDLT, 3, 3> LDLT3f;

/** LDLT, 4x4 float matrix decomposition method short definition */
typedef Decompose<float, DecompositionMethod::LDLT, 4, 4> LDLT4f;

/** LDLT, 3x3 double matrix decomposition method short definition */
typedef Decompose<double, DecompositionMethod::LDLT, 3, 3> LDLT3d;

/** LDLT, 4x4 double matrix decomposition method short definition */
typedef Decompose<double, DecompositionMethod::LDLT, 4, 4> LDLT4d;

/** LDLT decomposition */
template <typename S, int R, int C> class Decompose<S, DecompositionMethod::LDLT, R, C> {
public:

    typedef internal::SMatrix<S, R, C> Matrix;
    typedef internal::SMatrix<S, R, 1> Vector;
    typedef internal::Operator<S, internal::Operation<S, Matrix, Matrix, internal::OperationsType::M_TRANSPOSE>> Transpose;

    /** constructor */
    Decompose(const internal::SMatrix<S, R, C>& m) : M(m), L(std::string("L")), D(std::string("D")) {
        static_assert(internal::traits<S>::rows==1 && internal::traits<S>::cols==1, "Dimension error");
        static_assert(R==C and R>=3, "Dimension error");
        decompose();
    }

    /** Get the L matrix */
    inline const internal::SMatrix<S, R, C>& getL() const {
        return L;
    }

    /** Get the D matrix */
    inline const internal::SMatrix<S, R, C>& getD() const {
        return D;
    }

    /** Get the P matrix */
    inline internal::SMatrix<S, R, R> getPMatrix() const {
        /* get the P matrix from the P vector */
        internal::SMatrix<S,R,C> ret = internal::SMatrix<S,R,C>::identity();
        unsigned int r0=0, r1=0;
        for(unsigned int k=0; k<R-1; k++){
            r0 = k; P[k];
            ret.switchRows(r0<r1?r0:r1, r0<r1?r1:r0);
        }

        return ret;
    }

    /** solve the systema using the LDLT decomposition */
    inline internal::SMatrix<S, R, 1> solve (const internal::SMatrix<S, R, 1>& b) const {
        assert(decomposed);

        // Solve Ax = b given P A Pt = L D Lt
        //
        // L w = P b
        // D y = w
        // Lt z = y
        // x = Pt z
        //

        Matrix Per = getPMatrix();
        Vector Pb = Per * b;
        Vector w;
        Solve<S, SolveMethod::FORWARD>::template solve<Matrix, Vector, Vector>(L, w, Pb);
        Vector y;
        Solve<S, SolveMethod::DIRECT>::template solve<Matrix, Vector, Vector>(D, y, w);
        Vector z;
        Solve<S, SolveMethod::BACK>::template solve<Transpose, Vector, Vector>(L.transpose(), z, y);
        Vector x = Per.transpose() * z;
        return x;
    }

    /** Get P vector */
    inline const internal::SMatrix<S, R, 1>& getPVector() const {
        return P;
    }

    /** Check if the matrix is decomposed */
    inline bool isDecomposed() const {
        return decomposed;
    }

protected:

    inline void decompose() {
        decomposed = _LDLT_Implementation<S, _LDLT_Method_Type::PARTIAL_PIVOTING>::template compute<internal::SMatrix<S,R,C>,R,C>(M, P, L, D);
    }

    const internal::SMatrix<S, R, C>& M;
    internal::SMatrix<S, R, C> L;
    internal::SMatrix<S, R, C> D;
    internal::SMatrix<S, R, 1> P;
    bool decomposed = false;
};


} // solve

} // bsnlib

#endif // LDLT_H
