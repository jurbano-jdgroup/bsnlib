#ifndef SARRAYOBJECT_H
#define SARRAYOBJECT_H

#include <stdio.h>
#include <assert.h>

namespace bsnlib {

namespace internal {

// TODO: do align macros
// TODO: gcc and linux macros

/** static array object */
template <typename Scalar, int Dims, int Rows, int Cols> class SArrayObject {
public:

    /** default constructor */
    SArrayObject() {
        static_assert(traits<Scalar>::rows==1 && traits<Scalar>::cols==1, "Static array, scalar dimension error");
    }

    /** other element based constructor */
    SArrayObject(const SArrayObject& other){
        if(this != &other){
            for(unsigned int i=0; i<(Dims*Cols*Rows); ++i){
                (*this)[i] = other[i];
            }
        }
    }

    /** index operator */
    inline Scalar& operator[](const unsigned int index){
        assert(index < (Dims * Cols * Rows));
        return __data[index];
    }

    /** const index operator */
    inline const Scalar& operator[](const unsigned int index) const {
        assert(index < (Dims * Cols * Rows));
        return __data[index];
    }

    /** entry access operator */
    inline Scalar& operator()(const unsigned int dim, const unsigned int i, const unsigned int j){
        assert((dim<Dims) && (i<Rows) && (j<Cols));
        return __data[(dim*(Cols*Rows)) + (i*Cols) + j];
    }

    /** const entry access operator */
    inline const Scalar& operator()(const unsigned int dim, const unsigned int i, const unsigned int j) const {
        assert((dim<Dims) && (i<Rows) && (j<Cols));
        return __data[(dim*(Cols*Rows)) + (i*Cols) + j];
    }

    /** multiply each entry by a scalar value */
    template <typename Scalar2> inline void productByScalarEach(const Scalar2& val){
        static_assert(traits<Scalar2>::rows==1 && traits<Scalar2>::cols==1, "array object product by scalar, scalar2 type dimension error");
        for(unsigned int i=0; i<(Dims*Cols*Rows);++i){
            (*this)[i] *= static_cast<Scalar>(val);
        }
    }

    /** add a scalar value to each element */
    template<typename Scalar2> inline void addScalarEach(const Scalar2& val){
        static_assert(traits<Scalar2>::rows==1 && traits<Scalar2>::cols==1, "array object add scalar, scalar2 type dimension error");
        for(unsigned int i=0; i<(Dims*Cols*Rows);++i){
            (*this)[i] += static_cast<Scalar>(val);
        }
    }

    /** substract a scalar value to each element */
    template<typename Scalar2> inline void subsScalarEach(const Scalar2& val){
        static_assert(traits<Scalar2>::rows==1 && traits<Scalar2>::cols==1, "array object subs scalar, scalar2 type dimension error");
        for(unsigned int i=0; i<(Dims*Cols*Rows);++i){
            (*this)[i] -= static_cast<Scalar>(val);
        }
    }

    /** get the mex value of the entries */
    inline Scalar array_max() const {
        Scalar tmp = static_cast<Scalar>(0);
        for (int j=0; j<Dims * Rows * Cols; ++j) {
            if(tmp < (*this)[j]){
                tmp = (*this)[j];
            }
        }

        return tmp;
    }

    /** get the index of the max value of the entries */
    inline unsigned int array_max_index() const {
        Scalar tmp = static_cast<Scalar>(0);
        unsigned int index = 0;
        for (unsigned int j=0; j<Dims * Rows * Cols; ++j) {
            if(tmp < (*this)[j]){
                tmp = (*this)[j];
                index = j;
            }
        }

        return index;
    }

    /** get the data pointer */
    Scalar* data(){
        return __data;
    }

    /** get const data pointer */
    const Scalar* data() const {
        return __data;
    }

protected:

    Scalar __data[Dims * Rows * Cols];
};

} // internal

} // bsnlib

#endif // SARRAYOBJECT_H
