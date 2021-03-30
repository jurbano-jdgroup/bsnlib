
#ifndef TRAITS_H
#define TRAITS_H

namespace bsnlib {

namespace internal {

/** matrix traits */
template <class Matrix> struct traits {
    enum{
        rows = 1,
        cols = 1
    };
};

/** static matrix object definition */
template <typename Scalar, int Rows, int Cols> class SMatrixObject;

/** static matrix object traits */
template <typename Scalar, int Rows, int Cols> struct traits<SMatrixObject<Scalar, Rows, Cols>>{
    enum{
        rows = Rows,
        cols = Cols
    };
};

/** const static matrix object traits */
template <typename Scalar, int Rows, int Cols> struct traits<const SMatrixObject<Scalar, Rows, Cols>>{
    enum{
        rows = Rows,
        cols = Cols
    };
};

/** static matrix reference object traits */
template <typename Scalar, int Rows, int Cols> struct traits<SMatrixObject<Scalar, Rows, Cols>&>{
    enum{
        rows = Rows,
        cols = Cols
    };
};

/** const static matrix reference object traits */
template <typename Scalar, int Rows, int Cols> struct traits<const SMatrixObject<Scalar, Rows, Cols>&>{
    enum{
        rows = Rows,
        cols = Cols
    };
};

/** is scalar trait function */
template <typename _type> struct is_scalar {
    static constexpr bool value = false;
};

#define BSN_MAKE_IS_SCALAR_FUN(_obj_type) \
    template <> struct is_scalar<_obj_type> { \
    static constexpr bool value = true; \
}

BSN_MAKE_IS_SCALAR_FUN(float);
BSN_MAKE_IS_SCALAR_FUN(const float);
BSN_MAKE_IS_SCALAR_FUN(float&);
BSN_MAKE_IS_SCALAR_FUN(const float&);
BSN_MAKE_IS_SCALAR_FUN(double);
BSN_MAKE_IS_SCALAR_FUN(const double);
BSN_MAKE_IS_SCALAR_FUN(double&);
BSN_MAKE_IS_SCALAR_FUN(const double&);

// ////////////////////////////////
// enable complex values
// ////////////////////////////////

/** general float */
template <> struct is_scalar<std::complex<float>> {
    static constexpr bool value = true;
};

/** const general float */
template <> struct is_scalar<const std::complex<float>> {
    static constexpr bool value = true;
};

/** ref float */
template <> struct is_scalar<std::complex<float>&> {
    static constexpr bool value = true;
};

/** const ref float */
template <> struct is_scalar<const std::complex<float>&> {
    static constexpr bool value = true;
};

/** general float */
template <> struct is_scalar<std::complex<double>> {
    static constexpr bool value = true;
};

/** const general float */
template <> struct is_scalar<const std::complex<double>> {
    static constexpr bool value = true;
};

/** ref float */
template <> struct is_scalar<std::complex<double>&> {
    static constexpr bool value = true;
};

/** const ref float */
template <> struct is_scalar<const std::complex<double>&> {
    static constexpr bool value = true;
};

} // internal

} // bsnlib

#endif // TRAITS_H
