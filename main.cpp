#include "bsnlib.h"

using namespace bsnlib::internal;
using namespace bsnlib::solve;

int main() {
    SMatrix4d mat = SMatrix4d::random();

    mat.dump();

    DSlice<double, SMatrix4d&> slice(mat, 0, 2, 0, 2);
    slice.dump();

    std::cout<<"rows: "<<DSlice<double, SMatrix4d&>::_BaseRows<<std::endl;
}