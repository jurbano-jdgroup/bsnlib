#ifndef BSNLIB_H
#define BSNLIB_H

// check c++ 11
#ifdef __cplusplus
    #if __cplusplus < 201103L
        #error "C++ 11 required"
    #endif
#else
    #error "C++ required"
#endif

// platform and arch checking
#if defined (_WIN32) || defined (_WIN64) || defined (__WIN32__) || defined (__WIN64__)
    #define BSN_WIN 1
    #if defined (_WIN64) || defined (__WIN64__) || defined (__amd64__) || defined (__amd64)
        #define BSN_64 1
    #else
        #define BSN_32 1
    #endif
    
#elif defined(__linux__) || defined(__ANDROID__)
    #define BSN_LINUX
    #if defined(__amd64__) || defined(__aarch64__) || defined(__x86_64__)
        #define BSN_64 1
    #else
        #define BSN_32
    #endif
#endif

#if defined (__GNUG__) || defined (__MINGW32__) || defined (__MINGW64__)
#if !defined (BSN_64) && (defined (__x86_64__) || defined (__aarch64__))
#define BSN_64 1
#endif
#endif

// default 32 bits system
#if !defined (BSN_64) && !defined (BSN_32)
#define BSN_32
#endif

// check c++11 fatures
//#if defined (__cpp_constexpr) && defined (__cpp_decltype) && defined (__cpp_static_assert) && defined (__cpp_exceptions)
//#else
//#error "some features required"
//#endif

#define BSNLIB_DEBUG = 0
#define BSN_MAKE_DEBUG(msg) {if(BSNLIB_DEBUG){printf("%s", msg);}}

#define BSNLIB_EPS = 1e-6

// include util headers
#include <string>
#include <stdio.h>
#include <assert.h>
#include <limits>
#include <complex>
#include <math.h>

// include core headers
#include "bsnlib/internal/traits.h"
#include "bsnlib/internal/complex.h"
#include "bsnlib/internal/sarrayobject.h"
#include "bsnlib/internal/smatrixobject.h"
#include "bsnlib/internal/slice.h"
#include "bsnlib/internal/dslice.h"
#include "bsnlib/internal/dot.h"
#include "bsnlib/internal/raw.h"
#include "bsnlib/internal/operator.h"
#include "bsnlib/internal/add.h"
#include "bsnlib/internal/subs.h"
#include "bsnlib/internal/multiply.h"
#include "bsnlib/internal/smultiply.h"
#include "bsnlib/internal/transpose.h"
#include "bsnlib/internal/determinant.h"
#include "bsnlib/internal/decompose.h"
#include "bsnlib/internal/basern.h"
#include "bsnlib/internal/smatrix.h"

#include "bsnlib/solve/solve.h"
#include "bsnlib/solve/backsubstitution.h"
#include "bsnlib/solve/forwardsubstitution.h"
#include "bsnlib/solve/directsubstitution.h"
#include "bsnlib/solve/plu.h"
#include "bsnlib/solve/ldlt.h"
#include "bsnlib/solve/householder.h"
#include "bsnlib/solve/givens.h"
#include "bsnlib/solve/qr.h"
#include "bsnlib/solve/jacobi.h"
#include "bsnlib/solve/eigen.h"
#include "bsnlib/solve/svd.h"

#include "bsnlib/geometry/so2.h"
#include "bsnlib/geometry/so3.h"
#include "bsnlib/geometry/se3.h"

#endif // BSNLIB_H