
#ifndef SE3_H
#define SE3_H

namespace bsnlib {

namespace geometry {

/** SE3 definition */
template <typename S> class Se3;

/** Float SE3 short definition */
typedef Se3<float> Se3f;

/** Double SE3 short definition */
typedef Se3<double> Se3d;

/** SE3 implementation */
template <typename S> class Se3 {
public:

    /** default constructor */
    Se3() : so3(), t(std::string("t")) {
        assert(internal::traits<S>::rows==1 && internal::traits<S>::cols==1);
    }

    /** Construct this element based on another SE3 object */
    Se3(const Se3& other) {
        if(this!=&other) {
            so3 = other.getSo3();
            t = other.getT();
        }
    }

    /** Construct this element based on a so3 and t objects specified */
    Se3(const So3<S>& _so3, const internal::SMatrix<S,3,1>& _t) : so3(_so3), t(_t) {
        assert(internal::traits<S>::rows==1 && internal::traits<S>::cols==1);
    }

    /** Construct this element based on a rotation matrix and a traslation vector */
    Se3(const internal::SMatrix<S,3,3>& _rot, const internal::SMatrix<S,3,1>& _t): so3(_rot), t(_t) {
        assert(internal::traits<S>::rows==1 && internal::traits<S>::cols==1);
    }

    /** Equal copy operator */
    Se3& operator = (const Se3& other) {
        if(this!=&other) {
            so3 = other.getSo3();
            t = other.getT();
        }

        return (*this);
    }

    /** Get the so3 */
    inline So3<S>& getSo3() {
        return so3;
    }

    /** Get the so3*/
    inline const So3<S>& getSo3() const {
        return so3;
    }

    /** Get the traslation vector */
    inline const internal::SMatrix<S, 3, 1>& getT() const {
        return t;
    }

    /** Get the traslation vector */
    inline internal::SMatrix<S,3,1>& getT() {
        return t;
    }

    /** Get the inverse of this element */
    inline Se3 inverse() const {
        return Se3<S>(so3.getR().transpose(), so3.getR()*t*-1);
    }

    /** Multiply this element by a vector */
    inline internal::SMatrix<S,3,1> operator * (const internal::SMatrix<S,3,1>& x) const {
        return internal::SMatrix<S,3,1>(so3.getR()*x + t);
    }

    /** Concatenate */
    inline Se3& operator *= (const Se3<S>& other) {
        t += so3 * other.getT();
        so3 *= other.getSo3();

        return (*this);
    }

    /** Multiply this element by another So3 element */
    inline Se3 operator * (const Se3<S>& other) const {
        return Se3<S>(so3*other.getSo3(), t + (so3*other.getT()));
    }

    /** Get the skew vector representation of the SE3 */
    static internal::SMatrix<S,6,1> ln(const Se3& _se3) {
        internal::SMatrix<S,6,1> ret;

        const internal::SMatrix<S,3,1> _w = So3<S>::ln(_se3.getSo3());
        const S theta = std::sqrt(_w.dot(_w));

        ret.template slice<3,5,0,0>() = _w;
        internal::SMatrix<S,3,3> w_hat = So3<S>::hat(_w);

        if(theta<1e-10) {
            internal::SMatrix<S,3,3> g_inv = internal::SMatrix<S,3,3>::identity() - (w_hat*0.5) +
                    (w_hat*w_hat*(1./12.));

            ret.template slice<0,2,0,0>() = g_inv * _se3.getT();
        }
        else  {
            internal::SMatrix<S,3,3> g_inv = internal::SMatrix<S,3,3>::identity() - (w_hat*0.5) +
                    (w_hat*w_hat*((1.-(theta/(2*std::tan(theta/2.))))/(theta*theta)));

            ret.template slice<0,2,0,0>() = g_inv * _se3.getT();
        }

        return ret;
    }

    /** Get the SE3 from the skew vector u */
    static Se3 exp(const internal::SMatrix<S,6,1>& _u) {
        internal::SMatrix<S,3,1> u;
        u[0]=_u[0]; u[1]=_u[1]; u[2]=_u[2];

        internal::SMatrix<S,3,1> w;
        w[0]=_u[3]; w[1]=_u[4]; w[2]=_u[5];

        const S theta = std::sqrt(w.dot(w));

        const internal::SMatrix<S,3,3> w_hat = So3<S>::hat(w);

        internal::SMatrix<S,3,1> t;

        So3<S> _so3 = So3<S>::exp(w);

        if(theta<1e-10) {
            t = _so3*u;
        }
        else {
            internal::SMatrix<S,3,3> g;
            g = internal::SMatrix3d::identity() + (w_hat*((1-std::cos(theta))/(theta*theta))) +
                    (w_hat*w_hat*((theta-std::sin(theta))/(theta*theta*theta)));

            t = g*u;
        }

        return Se3<S>(_so3, t);
    }

protected:

    /** so3 */
    So3<S> so3;

    /** traslation vector */
    internal::SMatrix<S, 3, 1> t;
};

} // geometry

} // bsnlib

#endif // SE3_H
