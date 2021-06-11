
#ifndef D_SLICE_H
#define D_SLICE_H

//#include "bsnlib/internal/traits.h"
#include <type_traits>

namespace bsnlib {

namespace internal {

/** Slice entry update struct */
template <typename S, typename _MatrixB, bool _scalarType> struct DSliceEntryUpdate {
    static void update(S& A, const _MatrixB& value, unsigned int i, unsigned int j) {
        A = value;
    }
};

/** false specialization */
template <typename S, typename _MatrixB> struct DSliceEntryUpdate<S, _MatrixB, false> {
    static void update(S& A, const _MatrixB& value, unsigned int i, unsigned int j) {
        A = value(i, j);
    }
};

/** Dynamic slice definition */
template <typename S, typename _Matrix> class DSlice {
public:
    typedef _Matrix _BaseMatrix;
    typedef S _BaseType;
    static constexpr unsigned int _BaseRows = traits<_Matrix>::rows;
    static constexpr unsigned int _BaseCols = traits<_Matrix>::cols;

    /** Constructor */
    DSlice(_Matrix other, unsigned int RI, unsigned int RS,
           unsigned int CI, unsigned int CS): RI(RI), RS(RS), CI(CI), CS(CS), matrix(other) {
        static_assert(!std::is_pointer<_Matrix>::value, "Pointers not allowed here");
        static_assert(traits<S>::rows==1 && traits<S>::cols==1, "Dimension error");

        // update rows and cols
        rows = static_cast<unsigned int>(traits<_Matrix>::rows);
        cols = static_cast<unsigned int>(traits<_Matrix>::cols);

        // check RI, ...
        assert(RS<rows && RI<=RS);
        assert(CS<cols && CI<=CS);
    }

    /** Non constant entry access operator */
    S& operator()(unsigned int i, unsigned int j) {
        static_assert (!std::is_const<typename std::remove_reference<_Matrix>::type>::value,
                       "Mutable access for constant reference");
        assert(i<=(RS-RI) && j<=(CS-CI));
        return matrix(RI+i, CI+j);
    }

    /** Constant entry access operator */
    const S& operator()(unsigned int i, unsigned int j) const  {
        assert(i<=(RS-RI) && j<=(CS-CI));
        return matrix(RI+i, CI+j);
    }

    /** Equal operator */
    template <typename _Matrix2> DSlice& operator = (const _Matrix2& other) {
        static_assert (!std::is_const<typename std::remove_reference<_Matrix>::type>::value,
                       "Mutable access for constant reference");
        // dimensions
        constexpr bool isScalar = is_scalar<_Matrix2>::value;
        assert ((traits<_Matrix2>::rows==(RS-RI+1) && traits<_Matrix2>::cols==(CS-CI+1)) || isScalar);
        for(unsigned int i=0; i<(RS-RI+1); ++i){
            for(unsigned int j=0; j<(CS-CI+1); ++j){
                SliceEntryUpdate<S, _Matrix2, isScalar>::update(matrix(RI+i,CI+j), other, i, j);
            }
        }
        return (*this);
    }

    /** Display the name and entries of this matrix object */
    void dump() const {
        std::string tmpName = "matrix object";

        printf("\n[%s]:\n", tmpName.c_str());
        for (unsigned int i=0;i<(RS-RI+1);++i) {
            printf("[");
            for (unsigned int j=0;j<(CS-CI+1);++j) {
                printf("%.4f\t", (*this)(i,j));
            }
            printf("]\n");
        }
    }

protected:
    unsigned int RI, RS, CI, CS;
    unsigned int rows;
    unsigned int cols;
    _Matrix matrix;
};

} // internal

} // bsnlib

#endif // D_SLICE_H
