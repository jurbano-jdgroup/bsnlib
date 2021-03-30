#ifndef SMATRIX_H
#define SMATRIX_H

#include <math.h>
#include <complex>
#include <time.h>

/*#include "bsnlib/internal/traits.h"
#include "bsnlib/internal/complex.h"
#include "bsnlib/internal/slice.h"
#include "bsnlib/internal/dot.h"
#include "bsnlib/internal/operator.h"
#include "bsnlib/internal/sarrayobject.h"
#include "bsnlib/internal/smatrixobject.h"
#include "bsnlib/internal/transpose.h"
#include "bsnlib/internal/decompose.h"
#include "bsnlib/internal/determinant.h"*/

namespace bsnlib {

namespace internal {

/** SMatrix type definition */
template <typename _S, int _R, int _C> class SMatrix;

/** Static matrix traits definition */
template <typename Scalar, int Rows, int Cols> struct traits<SMatrix<Scalar, Rows, Cols>> {
    enum {
        rows = Rows,
        cols = Cols
    };
};

/** Static const matrix traits definition */
template <typename Scalar, int Rows, int Cols> struct traits<const SMatrix<Scalar, Rows, Cols>> {
    enum {
        rows = Rows,
        cols = Cols
    };
};

/** Static matrix ref traits definition */
template <typename Scalar, int Rows, int Cols> struct traits<SMatrix<Scalar, Rows, Cols>&> {
    enum {
        rows = Rows,
        cols = Cols
    };
};

/** Static matrix const ref traits definition */
template <typename Scalar, int Rows, int Cols> struct traits<const SMatrix<Scalar, Rows, Cols>&> {
    enum {
        rows = Rows,
        cols = Cols
    };
};

/**
 * @brief SMatrix2d
 */
typedef SMatrix<double, 2, 2> SMatrix2d;

/**
 * @brief SMatrix3d
 */
typedef SMatrix<double, 3, 3> SMatrix3d;

/**
 * @brief SMatrix4d
 */
typedef SMatrix<double, 4, 4> SMatrix4d;

/**
 * @brief SMatrix2f
 */
typedef SMatrix<float, 2, 2> SMatrix2f;

/**
 * @brief SMatrix3f
 */
typedef SMatrix<float, 3, 3> SMatrix3f;

/**
 * @brief SVector2d
 */
typedef SMatrix<double, 2, 1> SVector2d;

/**
 * @brief SVector3d
 */
typedef SMatrix<double, 3, 1> SVector3d;

/**
 * @brief SVector2f
 */
typedef SMatrix<float, 2, 1> SVector2f;

/**
 * @brief SVector3f
 */
typedef SMatrix<float, 3, 1> SVector3f;

/** Implementation  of the static matrix object */
template <typename Scalar, int Rows, int Cols> class SMatrix : public SMatrixObject<Scalar, Rows, Cols>, public SArrayObject<Scalar, 1, Rows, Cols> {
public:

    typedef SMatrix<Scalar, Rows, Cols> __type;
    typedef const SMatrix<Scalar, Rows, Cols> __const_type;
    typedef __type& __type_ref;
    typedef const __type& __const_type_ref;

    typedef SMatrixObject<Scalar, Rows, Cols> _BaseMatrix;
    typedef SArrayObject<Scalar, 1, Rows, Cols> Base;
    using Base::operator();
    using Base::data;

    ////////////////////////////////////////////////////
    /// Construction functions
    ////////////////////////////////////////////////////

    /** Default constructor */
    SMatrix(): _BaseMatrix(), Base() {}

    /** Constructor based on a string name */
    SMatrix(const std::string& _name) : _BaseMatrix(_name), Base() {}

    /** constructor based on another matrix object */
    SMatrix(const SMatrix<Scalar, Rows, Cols>& other) : _BaseMatrix(), Base() {
        if(this != &other){
            for(int i=0; i<Rows; ++i){
                for(int j=0; j<Cols; ++j){
                    (*this)(i,j) = other(i,j);
                }
            }
        }
    }

    /** Constructor based on another base matrix object */
    SMatrix(const SMatrixObject<Scalar, Rows, Cols>& other) : _BaseMatrix(), Base() {
        if(this != &other){
            for(int i=0; i<Rows; ++i){
                for(int j=0; j<Cols; ++j){
                    (*this)(i,j) = other(i,j);
                }
            }
        }
    }

    /** Operator based constructor */
    template <typename _A, typename _B, int _Op> SMatrix(const Operator<Scalar, Operation<Scalar, _A, _B, _Op>>& other) : _BaseMatrix(), Base() {
        static_assert(int(traits<Operator<Scalar, Operation<Scalar, _A, _B, _Op>>>::rows)==Rows, "Dimension error");
        static_assert(int(traits<Operator<Scalar, Operation<Scalar, _A, _B, _Op>>>::cols)==Cols, "Dimension error");
        for(int i=0; i<Rows; ++i){
            for(int j=0; j<Cols; ++j){
                (*this)(i,j) = other(i,j);
            }
        }
    }

    /** Build a vector of 3 entries */
    SMatrix(const Scalar& v0, const Scalar& v1, const Scalar& v2) : _BaseMatrix(), Base() {
        /* check if a row or column vector */
        static_assert((Rows==3 && Cols==1) || (Cols==3 && Rows==1), "Dimension error");

        (*this).data()[0] = v0;
        (*this).data()[1] = v1;
        (*this).data()[2] = v2;
    }

    /** build a vector of two entries */
    SMatrix(const Scalar& v0, const Scalar& v1) : _BaseMatrix(), Base() {
        /* check if a row or column vector */
        static_assert((Rows==2 && Cols==1) || (Cols==2 && Rows==1), "Dimension error");

        (*this).data()[0] = v0;
        (*this).data()[1] = v1;
    }

    /** Get a identity matrix */
    static __type identity() {
        SMatrix<Scalar, Rows, Cols> ret("Identity");
        for (unsigned int i=0; i<Rows; ++i){
            for (unsigned int j=0; j<Cols; ++j){
                ret(i, j) = i==j ? static_cast<Scalar>(1) : static_cast<Scalar>(0);
            }
        }

        return ret;
    }

    ////////////////////////////////////////////////////
    /// Assignment constructor operations
    ////////////////////////////////////////////////////

    /** Equal assignment operator */
    template <typename _MatrixType> __type_ref operator = (const _MatrixType& other) {
        static_assert(traits<_MatrixType>::rows==Rows && traits<_MatrixType>::cols==Cols, "Equal operator. Dimension mismatch");
        for(int i=0; i<Rows; ++i){
            for(int j=0; j<Cols; ++j){
                (*this)(i,j) = other(i,j);
            }
        }
        return (*this);
    }

    ////////////////////////////////////////////////////
    /// Entries selection
    ////////////////////////////////////////////////////

    /** Entry selection operator */
    Scalar& operator()(const unsigned int i, const unsigned int j)     {
        assert(i<Rows && j<Cols);
        return Base::operator()(0,i,j);
    }

    /** Const entry selection operator */
    const Scalar& operator()(const unsigned int i, const unsigned int j) const     {
        assert(i<Rows && j<Cols);
        return Base::operator()(0,i,j);
    }

    /** Get the x component if vector */
    const Scalar& x() const {
        /* check the dimension */
        static_assert((Rows>1 && Cols==1) || (Cols>1 && Rows==1), "Dimension error");
        return (*this).data()[0];
    }

    /** Get the x component if vector */
    Scalar& x() {
        /* check the dimension */
        static_assert((Rows>1 && Cols==1) || (Cols>1 && Rows==1), "Dimension error");
        return (*this).data()[0];
    }

    /** Get the y component if vector */
    const Scalar& y() const {
        /* check the dimension */
        static_assert((Rows>=2 && Cols==1) || (Cols>=2 && Rows==1), "Dimension error");
        return (*this).data()[1];
    }

    /** Get the y component if vector */
    Scalar& y() {
        /* check the dimension */
        static_assert((Rows>=2 && Cols==1) || (Cols>=2 && Rows==1), "Dimension error");
        return (*this).data()[1];
    }

    /** Get the y component if vector */
    const Scalar& z() const {
        /* check the dimension */
        static_assert((Rows>=3 && Cols==1) || (Cols>=3 && Rows==1), "Dimension error");
        return (*this).data()[2];
    }

    /** Get the y component if vector */
    Scalar& z() {
        /* check the dimension */
        static_assert((Rows>=3 && Cols==1) || (Cols>=3 && Rows==1), "Dimension error");
        return (*this).data()[2];
    }

    /** Get the y component if vector */
    const Scalar& w() const  {
        /* check the dimension */
        static_assert((Rows>=4 && Cols==1) || (Cols>=4 && Rows==1), "Dimension error");
        return (*this).data()[3];
    }

    /** Get the y component if vector */
    Scalar& w() {
        /* check the dimension */
        static_assert((Rows>=4 && Cols==1) || (Cols>=4 && Rows==1), "Dimension error");
        return (*this).data()[3];
    }

    ////////////////////////////////////////////////////
    /// Matrix and scalar operations
    ////////////////////////////////////////////////////

    /** Addition operator */
    template <typename B> Operator<Scalar, Operation<Scalar, SMatrix, B, OperationsType::M_SUM_M>> operator + (const B& other) const {
        static_assert(traits<B>::rows == Rows && traits<B>::cols==Cols, "Dimension error");
        return Operator<Scalar, Operation<Scalar, SMatrix, B, OperationsType::M_SUM_M>>((*this), other);
    }

    /** Substraction operator */
    template <typename B>  Operator<Scalar, Operation<Scalar, SMatrix, B, OperationsType::M_SUB_M>> operator - (const B& other) const {
        static_assert(traits<B>::rows == Rows && traits<B>::cols==Cols, "Dimension error");
        return Operator<Scalar, Operation<Scalar, SMatrix, B, OperationsType::M_SUB_M>>((*this), other);
    }

    /** Multiplication operator */
    template <typename B>
    Operator<Scalar, Operation<Scalar, SMatrix, B,
    traits<B>::rows==1?(traits<B>::cols==1?OperationsType::M_MUL_S:
                       ((Cols==1&&traits<B>::cols==Rows)?OperationsType::M_MUL_M:-1)):
                       (traits<B>::rows==Cols?OperationsType::M_MUL_M:-1)>>
    operator * (const B& other) const
    {
        return Operator<Scalar, Operation<Scalar, SMatrix, B,
         traits<B>::rows==1?(traits<B>::cols==1?OperationsType::M_MUL_S:
         ((Cols==1&&traits<B>::cols==Rows)?OperationsType::M_MUL_M:-1)):
         (traits<B>::rows==Cols?OperationsType::M_MUL_M:-1)>>((*this), other);
    }

    /** Self addition operator */
    template <typename B> SMatrix& operator += (const B& other) {
        static_assert(traits<B>::rows==Rows && traits<B>::cols==Cols, "Dimension error");
        for (int i=0; i<Rows; ++i){
            for (int j=0; j<Cols; ++j){
                (*this)(i, j) += other(i, j);
            }
        }

        return (*this);
    }

    /** Self substraction operator */
    template <typename B> SMatrix& operator -= (const B& other) {
        static_assert(traits<B>::rows==Rows && traits<B>::cols==Cols, "Dimension error");
        for (int i=0; i<Rows; ++i){
            for (int j=0; j<Cols; ++j){
                (*this)(i, j) -= other(i, j);
            }
        }

        return (*this);
    }

    /** Self multiplication operator */
    template <typename B> SMatrix& operator *= (const B& other)    {
        static_assert((traits<B>::rows==Rows && traits<B>::cols==Cols) || (traits<B>::rows==1 && traits<B>::cols==1), "Dimension error");
        constexpr const bool isScalar = traits<B>::rows==1 && traits<B>::cols==1;
        const SMatrix tmp = (*this);
        (*this) = Operator<Scalar, Operation<Scalar, SMatrix, B,
                isScalar?(OperationsType::M_MUL_S):(OperationsType::M_MUL_M)>>(tmp, other);

        return (*this);
    }

    /** Transpose operation */
    Operator<Scalar, Operation<Scalar, SMatrix, SMatrix, OperationsType::M_TRANSPOSE>> transpose() const {
        return Operator<Scalar, Operation<Scalar, SMatrix, SMatrix, OperationsType::M_TRANSPOSE>>((*this));
    }

    ////////////////////////////////////////////////////
    /// Decompositions
    ////////////////////////////////////////////////////

    /** PLU decomposition */
    solve::Decompose<Scalar, solve::DecompositionMethod::PLU, Rows, Cols> plu() const {
        return solve::Decompose<Scalar, solve::DecompositionMethod::PLU, Rows, Cols>(*this);
    }

    /** LDLT decomposition object */
    solve::Decompose<Scalar, solve::DecompositionMethod::LDLT, Rows, Cols> ldlt() const {
        return solve::Decompose<Scalar, solve::DecompositionMethod::LDLT, Rows, Cols>(*this);
    }

    /** QR decomposition object */
    solve::Decompose<Scalar, solve::DecompositionMethod::QR, Rows, Cols> qr() const {
        return solve::Decompose<Scalar, solve::DecompositionMethod::QR, Rows, Cols>((*this), true);
    }

    // TODO: SVD, Eigen

    ////////////////////////////////////////////////////
    /// slice
    ////////////////////////////////////////////////////

    /** SMatrix slice */
    template <int RI, int RS, int CI, int CS> Slice<Scalar, __type&, RI, RS, CI, CS> slice() {
        return Slice<Scalar, __type&, RI, RS, CI, CS>((*this));
    }

    /** SMatrix slice */
    template <int RI, int RS, int CI, int CS> Slice<Scalar, __const_type&, RI, RS, CI, CS> slice() const {
        return Slice<Scalar, __const_type&, RI, RS, CI, CS>((*this));
    }

    ////////////////////////////////////////////////////
    /// Misc
    ////////////////////////////////////////////////////

    /** Set all the matrix entries to zero */
    void setZero() {
        for(int i=0; i<Rows * Cols; ++i){
            (*this)[i] = static_cast<Scalar>(0);
        }
    }

    /** clear the entries of the mat matrix except diagonal */
    static void makeIdentity(SMatrix& mat) {
        for(unsigned int i=0; i<Rows; ++i){
            for(unsigned int j=0; j<Cols; ++j){
                mat(i, j) = i==j ? static_cast<Scalar>(1) : static_cast<Scalar>(0);
            }
        }
    }

    /** Calculate the norm if vector */
    double norm() const {
        /** check if colum vector */
        static_assert((Cols==1&&Rows!=1) || (Rows==1&&Cols!=1), "Dimension error");
        double tmp = 0.0;
        for(unsigned int i=0; i<Rows; ++i){
            for(unsigned int j=0; j<Cols; ++j){
                tmp += complex_abs<Scalar>::compute(complex_conjugate<Scalar>::compute((*this)(i,j)) * (*this)(i,j));
            }
        }

        return std::sqrt(tmp);
    }

    /** Dot product for vectors only */
    template <typename _VectorType> Scalar dot (const _VectorType& other) const {
        /** check if colum vector */
        static_assert((Cols==1&&Rows>1) && (traits<_VectorType>::cols==1 && traits<_VectorType>::rows>1), "Dimension error");
        return Dot<Scalar>::template compute<SMatrix, SMatrix>((*this), other);
    }

    /** Cross product for 3 entries vectors only */
    template <typename _VectorType> SMatrix<Scalar,Rows,Cols> cross (const _VectorType& other) const {
        /** check if column vector or row vector */
        static_assert((Cols==1&&Rows>1) && (traits<_VectorType>::cols==1 && traits<_VectorType>::rows>1), "Dimension error");
        SMatrix<Scalar,Rows,Cols> ret;

        ret(0,0) = (*this)(1,0)*other(2,0) - (*this)(2,0)*other(1,0);
        ret(1,0) = (*this)(2,0)*other(0,0) - (*this)(0,0)*other(2,0);
        ret(2,0) = (*this)(0,0)*other(1,0) - (*this)(1,0)*other(0,0);

        return ret;
    }

    /** Normalize only if vector */
    void normalize() {
        /** check if colum vector or row vector */
        static_assert((Cols==1&&Rows!=1) || (Rows==1&&Cols!=1), "Dimension error");
        const double norm = (*this).norm();
        assert(norm!=static_cast<Scalar>(0));
        for(unsigned int i=0; i<Rows; ++i){
            for(unsigned int j=0; j<Cols; ++j){
                (*this)(i, j) *= 1.0 / norm;
            }
        }
    }

    /** Switch the R0 and R1 rows */
    void switchRows(const unsigned int R0, const unsigned int R1)     {
        assert(R0<Rows && R1<Rows && R0<=R1);
        if(R0==R1){return;}

        Scalar tmp = static_cast<Scalar>(0);

        for (int j=0; j<Cols; ++j){
            tmp = (*this)(R0, j);
            (*this)(R0, j) = (*this)(R1, j);
            (*this)(R1, j) = tmp;
        }
    }

    /** Switch the Col0 with the Col1 */
    void switchCols(const unsigned int Col0, const unsigned int Col1){
        assert(Col0<Cols && Col1<Cols && Col0<=Col1);

        if(Col1==Col0){return;}
        Scalar tmp = static_cast<Scalar>(0);

        for (int j=0; j<Rows; ++j) {
            tmp = (*this)(j, Col0);
            (*this)(j, Col0) = (*this)(j, Col1);
            (*this)(j, Col1) = tmp;
        }
    }

    /** Get a vector with the diagonal elements of the matrix */
    SMatrix<Scalar, Rows<Cols?Rows:Cols, 1> diagonal() const {
        static_assert(Cols>1 || Rows>1, "Dimension error");
        constexpr int min = Rows<Cols? Rows : Cols;
        SMatrix<Scalar, min, 1> vector;

        for (int i=0; i<min; ++i) {
            vector(i, 0) = (*this)(i, i);
        }

        return vector;
    }

    /** Get the ith col of the matrix */
    SMatrix<Scalar, Rows, 1> col(const unsigned int i) const {
        assert(i<Cols);
        SMatrix<Scalar, Rows, 1> vector;

        for (int j=0; j<Rows; ++j) {
            vector(j, 0) = (*this)(j, i);
        }

        return vector;
    }

    /** Get the ith row of the matrix */
    SMatrix<Scalar, 1, Cols> row(const unsigned int i) const {
        assert(i < Rows);
        SMatrix<Scalar, 1, Cols> vector;

        for (int j=0; j<Cols; ++j) {
            vector(0, j) = (*this)(i, j);
        }

        return vector;
    }

    /** get the ith row of the matrix */
    template <typename _MatrixType, typename OtherScalar, int otherCols> static SMatrix<OtherScalar,1,otherCols> row (_MatrixType &A, const unsigned int i) {
        // check other scalar
        static_assert (is_scalar<OtherScalar>::value, "Other scalar is not a scalar type");
        // check _MatrixTYpe
        constexpr unsigned int cols = (unsigned int) traits<_MatrixType>::cols;
        static_assert (cols>=2 && otherCols>1 && otherCols<=cols, "OtherCols error comparing with matrix type cols");
        // check i
        assert(i<int(traits<_MatrixType>::rows));

        SMatrix<OtherScalar,1,otherCols> ret;
        for(unsigned int j=0; j<otherCols; ++j){
            ret(0,j) = A(i,j);
        }

        return ret;
    }

    /** Get a matrix with the rows and cols specified */
    template <int R0, int R1, int C0, int C1> SMatrix<Scalar, R1 - R0 + 1, C1 - C0 + 1> subMatrix() const {
        static_assert(R0>=0 && R1>=0 && C0>=0 && C1>=0, "Dimension error");
        static_assert(R1<Rows && C1<Cols && R0<=R1 && C0<=C1, "Dimension error");
        SMatrix<Scalar, R1 - R0 + 1, C1 - C0 + 1> Matrix;

        for (int i=0; i<=(R1 - R0); ++i) {
            for (int j=0; j<=(C1 - C0); ++j) {
                Matrix(i, j) = (*this)(i + R0, j + C0);
            }
        }

        return Matrix;
    }

    /** Get a permutation matrix from a permutation vector */
    template <typename _VectorType> static __type permutationMatrix(const _VectorType& P) {
        static_assert (traits<_VectorType>::rows==Rows && traits<_VectorType>::cols==1, "Dimension error in permutation vector");
        // TODO: check if valid permutation for non square matrix
        static_assert (Rows==Cols, "Permutation matrix implemented for square matrices");
        __type ret = __type::identity();
        unsigned int r0=0, r1=0;
        for(unsigned int k=0; k<Rows-1; k++){
            r0 = k;
            r1 = P(k,0);
            ret.switchRows(r0<r1?r0:r1, r0<r1?r1:r0);
        }

        return ret;
    }

    /** get a random values matrix */
    static __type random() {
        __type ret;
        srand(time(NULL));

        for(unsigned int i=0; i<Rows; ++i) {
            for(unsigned int j=0; j<Cols; ++j) {
                ret(i, j) = static_cast<Scalar>(rand() % 100);
            }
        }

        return ret;
    }

    /** Compute the determinant of this matrix
    only if the dimension is lower than 4 and it's
    square*/
    Scalar determinant() const     {
        /* check the matrix dimension */
        static_assert(Rows==Cols && Rows<=3, "Dimension error");
        return Determinant<Scalar>::template calculate<SMatrix<Scalar, Rows, Cols>>(*this);
    }
};

} // bsnlib

} // bsnlib

#endif // SMATRIX_H
