
#include "bsnlib.h"

using namespace bsnlib::internal;
using namespace bsnlib::solve;

int main() 
{
    using SMatrix4d = SMatrix<double,4,4>;
    using SVector4dx = SMatrix<std::complex<double>,4,1>;
    SMatrix4d mat = SMatrix4d::random();

    mat.dump();

    SVector4dx eigenValues;
    const bool flag = EigenQrd::template compute<SMatrix4d, SVector4dx>(mat, eigenValues);

    if (!flag) {
        std::cout<<"Error encontrando los valores propios"<<std::endl;
        return -1;
    }

    eigenValues.dump();
}