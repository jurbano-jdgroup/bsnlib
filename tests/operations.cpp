#include "./../bsnlib.h"
#define TEST_EPS 1e-8

using namespace bsnlib::internal;

int main() {
    SMatrix3d mat = SMatrix3d::random();
    const SMatrix3d orig = mat * 1;
    
    const SMatrix3d substraction = mat - orig;
    assert(substraction.array_sum_sqr() < TEST_EPS);

    const double sum = orig.array_sum();
    const double sum_sqr = orig.array_sum_sqr();
    const double scalar_value = 3.05;

    const SMatrix3d scalarMul = orig * scalar_value;
    assert(std::abs(sum*scalar_value - scalarMul.array_sum()) < TEST_EPS);

    const SMatrix3d dob = mat + mat;
    assert(std::abs(sum*2 - dob.array_sum()) < TEST_EPS);
}