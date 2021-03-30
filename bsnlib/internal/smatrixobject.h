
#ifndef SMATRIXOBJECT_H
#define SMATRIXOBJECT_H

//#include "bsnlib/internal/traits.h"
//#include "bsnlib/internal/complex.h"
//#include <string>
#include <iostream>

namespace bsnlib {

namespace internal {

/** SMatrix object */
template<typename Scalar, int Rows, int Cols> class SMatrixObject{
public:

    /** default constructor */
    SMatrixObject() {
        static_assert(traits<Scalar>::rows==1 && traits<Scalar>::cols==1, "Dimension error");
    }

    /** name based constructor */
    SMatrixObject(const std::string& _name) : name(_name){
        static_assert(traits<Scalar>::rows==1 && traits<Scalar>::cols==1, "Dimension error");
    }

    //////////////////////////////////////////////////////
    /// Assignment operations
    //////////////////////////////////////////////////////

    /** entries assignment operator */
    virtual Scalar& operator()(const unsigned int i, const unsigned int j) = 0;

    /** entries assignment operator */
    virtual const Scalar& operator()(const unsigned int i, const unsigned int j) const = 0;

    //////////////////////////////////////////////////////
    /// Misc
    //////////////////////////////////////////////////////

    /** display the name and entries of this matrix object */
    void dump() const {
        std::string tmpName = "matrix object";
        if(name.compare(std::string(""))!=0){
            tmpName = name;
        }

        printf("\n[%s]:\n", tmpName.c_str());
        for (int i=0;i<Rows;++i) {
            printf("[");
            for (int j=0;j<Cols;++j) {
                std::cout<<(*this)(i,j)<<"\t";
            }
            printf("]\n");
        }
    }

protected:
    std::string name;
};

} // internal

} // bsnlib

#endif // SMATRIXOBJECT_H
