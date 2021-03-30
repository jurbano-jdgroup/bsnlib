
#ifndef SO2_H
#define SO2_H

#include <math.h>

namespace bsnlib {

namespace geometry {

/**
 * So2 class definition
 *
 */
template <typename Scalar> class So2;

/**
 * @brief so2d Double specialization of so2
 *
 */
typedef So2<double> So2d;

/**
 * @brief So2f Float specialization of so2
 *
 */
typedef So2<float> So2f;

// TODO: Only float and double scalar values

/**
 * SO2
 *
 */
template <typename Scalar> class So2
{
public:

    /**
     * @brief So2 Construct this element
     * based on an angle param
     *
     * @param _angle
     */
    So2(const Scalar& _angle) : angle(_angle)
    {
        assert (internal::traits<Scalar>::rows==1 && internal::traits<Scalar>::cols==1);

        const double _cos = std::cos(_angle);

        const double _sin = std::sin(_angle);

        R(0, 0) = static_cast<Scalar>(_cos);

        R(1, 1) = static_cast<Scalar>(_cos);

        R(0, 1) = static_cast<Scalar>(-1 * _sin);

        R(1, 0) = static_cast<Scalar>(_sin);
    }

    /**
     * @brief So2 Construct this element based on
     * a matrix
     *
     * @param _R
     */
    So2(const internal::SMatrixObject<Scalar, 2, 2>& _R) : R(_R) {
        assert (internal::traits<Scalar>::rows==1 && internal::traits<Scalar>::cols==1);

        assert(_R(0, 0)==_R(1, 1) && _R(0, 1)==(-1.0 * _R(1, 0)));

        const double a0_2 = std::pow(_R(0, 0), 2);

        const double a1_2 = std::pow(_R(1, 0), 2);

        const double a_err = std::abs((a0_2 + a1_2) - 1.0);

        assert(a_err <= 0.000001);

        // TODO: check this
        this->angle = std::atan2(_R(1, 0), _R(0, 0));
    }

    /**
     * @brief So2 Construct this element based
     * on another SO2 object
     *
     * @param other
     */
    So2(const So2& other) :
        angle(other.getAngle()), R(other.getMatrix())
    {
        assert (internal::traits<Scalar>::rows==1 && internal::traits<Scalar>::cols==1);
    }

    /**
     * @brief operator = Equal operator
     *
     * @param other
     * @return
     */
    So2& operator = (const So2& other)
    {
        if(this!=&other)
        {
            angle = other.getAngle();
            R = other.getMatrix();
        }

        return *this;
    }

    /**
     * @brief getAngle get the angle
     *
     * @return
     */
    const double& getAngle() const
    {
        return angle;
    }

    /**
     * @brief getMatrix Get the matrix
     *
     * @return
     */
    const internal::SMatrix<Scalar, 2, 2>& getMatrix() const
    {
        return R;
    }

    /**
     * @brief log Get the 1-dimension vector (scalar)
     * representation of the so2 specified object
     *
     * @param other
     * @return
     */
    static Scalar ln(const So2& other)
    {
        return other.getAngle();
    }

    /**
     * @brief exp Get the so2 from the angle specified
     * @param angle
     * @return
     */
    static So2 exp(const Scalar& angle)
    {
        return So2(angle);
    }

    /**
     * @brief operator * Apply this rotation
     * to the specified vector
     *
     * @param vec
     * @return
     */
    inline internal::SMatrix<Scalar, 2, 1> operator *(const internal::SMatrix<Scalar, 2, 1>& vec) const
    {
        return internal::SMatrix<Scalar, 2, 1>(R * vec);
    }


    /** concatenate to this element */
    inline So2& operator *= (const So2& other)
    {
        R *= other.getMatrix();

        return (*this);
    }


    /** Multiply this element by another So2 element */
    inline So2 operator * (const So2& other) const
    {
        return So2(R * other.getMatrix());
    }

    /**
     * @brief inverse ge the inverse of this element
     * @return
     */
    inline So2 inverse() const
    {
        return So2(R.transpose());
    }

protected:

    /**
     * @brief angle The angle of the rotation
     *
     */
    Scalar angle;

    /**
     * @brief R The matrix representation of the
     * rotation
     *
     */
    internal::SMatrix<Scalar, 2, 2> R;
};

} // geometry

} // bsnlib

#endif // SO2_H
