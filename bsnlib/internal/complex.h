#ifndef COMPLEX_H
#define COMPLEX_H

//#include <complex>

namespace bsnlib {

namespace internal {

/** is complex struct trait */
template <typename S> struct is_complex {
    static constexpr bool value =  false;
};

typedef std::complex<float> complexf;
typedef const std::complex<float> const_complexf;
typedef std::complex<float>& ref_complexf;
typedef const std::complex<float>& const_ref_complexf;

typedef std::complex<double> complexd;
typedef const std::complex<double> const_complexd;
typedef std::complex<double>& ref_complexd;
typedef const std::complex<double>& const_ref_complexd;

/** float complex */
template <> struct is_complex<std::complex<float>>{
    static constexpr bool value = true;
};

/** const float complex */
template <> struct is_complex<const std::complex<float>>{
    static constexpr bool value = true;
};

/** ref float complex */
template <> struct is_complex<std::complex<float>&>{
    static constexpr bool value = true;
};

/** const ref float complex */
template <> struct is_complex<const std::complex<float>&>{
    static constexpr bool value = true;
};

/** double complex */
template <> struct is_complex<std::complex<double>>{
    static constexpr bool value = true;
};

/** const double complex */
template <> struct is_complex<const std::complex<double>>{
    static constexpr bool value = true;
};

/** ref double complex */
template <> struct is_complex<std::complex<double>&>{
    static constexpr bool value = true;
};

/** const ref double complex */
template <> struct is_complex<const std::complex<double>&>{
    static constexpr bool value = true;
};

/** complex conjugate */
template <typename S> struct complex_conjugate {
    static S compute (const S& val) {
        return val;
    }
};

/** float complex conjugate */
template <> struct complex_conjugate<complexf>{
    static complexf compute (const complexf& val) {
        return std::conj(val);
    }
};

/** float complex conjugate */
template <> struct complex_conjugate<complexd>{
    static complexd compute (const complexd& val) {
        return std::conj(val);
    }
};

/** complex abs */
template <typename S> struct complex_abs {
    static S compute (const S& val) {
        return static_cast<S>(std::abs(val));
    }
};

/** float complex abs */
template <> struct complex_abs<complexf> {
    static float compute (const complexf& val){
        return std::abs(val);
    }
};

/** double complex abs */
template <> struct complex_abs<complexd> {
    static double compute (const complexd& val){
        return std::abs(val);
    }
};

} // internal

} // bsnlib

#endif // COMPLEX_H
