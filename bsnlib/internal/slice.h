
#ifndef SLICE_H
#define SLICE_H

//#include "bsnlib/internal/traits.h"
#include <type_traits>
#include "assert.h"
#include <string>

namespace bsnlib {

namespace internal {

/** Slice definition */
template <typename S, typename _MatrixType, int RI, int RS, int CI, int CS> class Slice;

/** Slice traits */
template <typename S, typename _MatrixType, int RI, int RS, int CI, int CS> struct traits<Slice<S, _MatrixType, RI, RS, CI, CS>> {
    enum{
        rows = RS - RI + 1,
        cols = CS - CI + 1
    };
};

/** Slice traits */
template <typename S, typename _MatrixType, int RI, int RS, int CI, int CS> struct traits<const Slice<S, _MatrixType, RI, RS, CI, CS>> {
    enum{
        rows = RS - RI + 1,
        cols = CS - CI + 1
    };
};

/** Slice traits */
template <typename S, typename _MatrixType, int RI, int RS, int CI, int CS> struct traits<Slice<S, _MatrixType, RI, RS, CI, CS>&> {
    enum{
        rows = RS - RI + 1,
        cols = CS - CI + 1
    };
};

/** Slice traits */
template <typename S, typename _MatrixType, int RI, int RS, int CI, int CS> struct traits<const Slice<S, _MatrixType, RI, RS, CI, CS>&> {
    enum{
        rows = RS - RI + 1,
        cols = CS - CI + 1
    };
};

/** Slice entry update struct */
template <typename S, typename _MatrixB, bool _scalarType> struct SliceEntryUpdate {
    static void update(S& A, const _MatrixB& value, unsigned int i, unsigned int j) {
        A = value;
    }
};

/** false specialization */
template <typename S, typename _MatrixB> struct SliceEntryUpdate<S, _MatrixB, false> {
    static void update(S& A, const _MatrixB& value, unsigned int i, unsigned int j) {
        A = value(i, j);
    }
};

/** Slice static implementation */
template <typename S, typename _MatrixType, int RI, int RS, int CI, int CS> class Slice {
public:

    /** constructor */
    Slice(_MatrixType other) : Matrix(other) {
        static_assert(!std::is_pointer<_MatrixType>::value, "No pointer allowed");
        static_assert(traits<S>::rows==1 && traits<S>::cols==1, "Dimension error");
        static_assert(RI<=RS && (RS <= traits<_MatrixType>::rows), "Dimension error");
        static_assert(CI<=CS && (CS <= traits<_MatrixType>::cols), "Dimension error");
    }

    /** other based constructor */
    //Slice(const Slice& other) {}

    /** entry acces operator */
    S& operator()(const unsigned int i, const unsigned int j) {
        static_assert (!std::is_const<typename std::remove_reference<_MatrixType>::type>::value, "Mutable access for constant reference");
        assert(i<=(RS-RI) && j<=(CS-CI));
        return Matrix(RI+i, CI+j);
    };

    /** entry acces operator */
    const S& operator()(const unsigned int i, const unsigned int j) const {
        assert(i<=(RS-RI) && j<=(CS-CI));
        return Matrix(RI+i, CI+j);
    };

    /** equal operator */
    template <typename _Matrix2> Slice& operator = (const _Matrix2& other) {
        static_assert (!std::is_const<typename std::remove_reference<_MatrixType>::type>::value, "Mutable access for constant reference");
        // dimensions
        constexpr bool isScalar = is_scalar<_Matrix2>::value;
        static_assert ((traits<_Matrix2>::rows==(RS-RI+1) && traits<_Matrix2>::cols==(CS-CI+1)) || isScalar, "Dimension error");
        for(unsigned int i=0; i<(RS-RI+1); ++i){
            for(unsigned int j=0; j<(CS-CI+1); ++j){
                SliceEntryUpdate<S, _Matrix2, isScalar>::update(Matrix(RI+i,CI+j), other, i, j);
            }
        }
        return (*this);
    }

    /** display the name and entries of this matrix object */
    void dump() const {
        std::string tmpName = "matrix object";

        printf("\n[%s]:\n", tmpName.c_str());
        for (int i=0;i<(RS-RI+1);++i) {
            printf("[");
            for (int j=0;j<(CS-CI+1);++j) {
                printf("%.4f\t", (*this)(i,j));
            }
            printf("]\n");
        }
    }

protected:

   /** Matrix pointer element */
   _MatrixType Matrix;
};

} // internal

} // bsnlib

#endif // SLICE_H
