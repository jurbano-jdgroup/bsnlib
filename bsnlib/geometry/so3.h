
#ifndef SO3_H
#define SO3_H

#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

namespace bsnlib {

namespace geometry {

/** SO3 definition */
template <typename S> class So3;

/**  Float so3 short definition */
typedef So3<float> So3f;

/** Double so3 short definition */
typedef So3<double> So3d;

/** Space rotation object */
template <typename S> class So3 {
public:

    /** Default constructor */
    So3() : R(internal::SMatrix<S, 3, 3>::identity()) {
        assert(internal::traits<S>::rows==1 && internal::traits<S>::cols==1);
    }

    /** Construct this element based on another so3 element */
    So3(const So3& other) {
        if(&other!=this) {
            const S det = other.getR().determinant();

            // TODO: assure the correct epsilon
            assert(det<static_cast<S>(1.001) && det>static_cast<S>(0.999));

            R = other.getR();
        }
    }

    /** Construct this element based on a rotation specified */
    So3(const internal::SMatrix<S,3,3>& _R) {
        const S det = _R.determinant();
        assert(det<static_cast<S>(1.001) && det>static_cast<S>(0.999));
        R = _R;
    }

    /** Equal operator */
    So3& operator = (const So3& other){
        if(&other!=this){
            const S det = other.getR().determinant();
            assert(det<static_cast<S>(1.001) && det>static_cast<S>(0.999));
            R = other.getR();
        }

        return (*this);
    }


    /** Get the rotation matrix */
    inline const internal::SMatrix<S, 3, 3>& getR() const {
        return R;
    }


    /** Multiply this element with the vector x specified */
    inline internal::SMatrix<S,3,1> operator * (const internal::SMatrix<S,3,1>& x) const {
        return internal::SMatrix<S,3,1>(R * x);
    }


    /** Concatenate to this element */
    inline So3& operator *= (const So3<S>& other) {
        const internal::SMatrix3d _this_R = R;
        R = other.getR() * _this_R;
        return (*this);
    }


    /** Multiply this element by another se3 element */
    inline So3 operator * (const So3<S>& other) const {
        return So3<S>(R * other.getR());
    }


    /** Get the inverse of this element */
    inline So3 inverse() const {
        return So3(R.transpose());
    }


    /** Get the so3 based on the rodrigues formula */
    static So3 exp(const internal::SMatrix<S,3,1>& w) {
        /* must use double precision */
        const double w_sqr = internal::Dot<S>::template compute<internal::SMatrix<S,3,1>, internal::SMatrix<S,3,1>>(w,w);
        const double tetha = std::sqrt(w_sqr);
        const double tetha_inv = 1.0 / tetha;

        const double one_6 = 1.0/6.0;
        const double one_20 = 1.0/20.0;

        double sin_t_t = 0.0;
        double cos_t_t2 = 0.0;

        if(w_sqr < 1e-8) {
            sin_t_t = 1.0 - one_6 * w_sqr;
            cos_t_t2 = 0.5;
        }
        else if(w_sqr < 1e-6) {
            cos_t_t2 = 0.5 - 0.25 * one_6 * w_sqr;
            sin_t_t = 1.0 - w_sqr * one_6*(1.0 - one_20 * w_sqr);
        }
        else{
            sin_t_t = std::sin(tetha) * tetha_inv;
            cos_t_t2 = (1.0 - std::cos(tetha)) * tetha_inv * tetha_inv;
        }

        internal::SMatrix3d ret;

        const double w0_2 = w[0] * w[0];
        const double w1_2 = w[1] * w[1];
        const double w2_2 = w[2] * w[2];

        ret(0, 0) = 1.0 - cos_t_t2*(w1_2 + w2_2);
        ret(0, 1) = (-1.0 * sin_t_t * w[2]) + (cos_t_t2 * w[1] * w[0]);
        ret(0, 2) = (w[1] * sin_t_t) + (cos_t_t2 * w[2] * w[0]);
        ret(1, 0) = (sin_t_t * w[2]) + (cos_t_t2 * w[0] * w[1]);
        ret(1, 1) = 1.0 - cos_t_t2*(w2_2 + w0_2);
        ret(1, 2) = (-1.0 * sin_t_t * w[0]) + (cos_t_t2 * w[2] * w[1]);
        ret(2, 0) = (-1.0 * sin_t_t * w[1]) + (cos_t_t2 * w[0] * w[2]);
        ret(2, 1) = (sin_t_t * w[0]) + (cos_t_t2 * w[2] * w[1]);
        ret(2, 2) = 1.0 - cos_t_t2*(w1_2 + w0_2);

        return So3(ret);
    }


    /** Get the exponential represention of an so3 element */
    static internal::SMatrix<S,3,1> ln(const So3& so3) {
        internal::SMatrix<S,3,1> result;

        const internal::SMatrix<S,3,3>& my_matrix = so3.getR();

        /* must use double precision */
        const double cos_angle = (my_matrix(0 ,0) + my_matrix(1, 1) + my_matrix(2, 2) - 1.0) * 0.5;
        result[0] = (my_matrix(2, 1) - my_matrix(1, 2))/2;
        result[1] = (my_matrix(0, 2) - my_matrix(2, 0))/2;
        result[2] = (my_matrix(1, 0) - my_matrix(0, 1))/2;

        S sin_angle_abs = std::sqrt(result.dot(result));

        if (cos_angle > M_SQRT1_2) {
            if(sin_angle_abs > 0){
                result.selfScalarMultiply(asin(sin_angle_abs) / sin_angle_abs);
            }
        }
        else if( cos_angle > - M_SQRT1_2) {
            const double angle = acos(cos_angle);
            result.selfScalarMultiply(angle / sin_angle_abs);
        }
        else {
            const double angle = M_PI - asin(sin_angle_abs);
            const double d0 = my_matrix(0, 0) - cos_angle,
                d1 = my_matrix(1, 1) - cos_angle,
                d2 = my_matrix(2, 2) - cos_angle;

            internal::SMatrix<S,3,1> r2;

            if(d0*d0 > d1*d1 && d0*d0 > d2*d2){
                r2[0] = d0;
                r2[1] = (my_matrix(1, 0)+my_matrix(0, 1))/2;
                r2[2] = (my_matrix(0, 2)+my_matrix(2, 0))/2;
            }
            else if(d1*d1 > d2*d2){
                r2[0] = (my_matrix(1, 0)+my_matrix(0, 1))/2;
                r2[1] = d1;
                r2[2] = (my_matrix(2, 1)+my_matrix(1, 2))/2;
            }
            else{
                r2[0] = (my_matrix(0, 2)+my_matrix(2, 0))/2;
                r2[1] = (my_matrix(2, 1)+my_matrix(1, 2))/2;
                r2[2] = d2;
            }

            if(r2.dot(result) < 0) {
                r2.selfScalarMultiply(-1.0);
            }

            r2.normalize();

            result = r2 * angle;
        }

        return result;
    }

    /** hat operator */
    static internal::SMatrix<S,3,3> hat(const internal::SMatrix<S,3,1>& _vec) {
        internal::SMatrix<S,3,3> ret;

        ret(0, 0) = 0; ret(0, 1) = _vec[2]*(-1); ret(0, 2) = _vec[1];
        ret(1, 0) = _vec[2]; ret(1, 1) = 0; ret(1, 2) = _vec[0]*(-1);
        ret(2, 0) = _vec[1]*(-1); ret(2, 1) = _vec[0]; ret(2, 2) = 0;

        return ret;
    }


protected:

    /** Rotation matrix */
    internal::SMatrix<S, 3, 3> R;
};

} // geometry

} // bsnlib

#endif // SO3_H
