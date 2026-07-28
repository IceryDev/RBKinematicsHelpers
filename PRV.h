#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>

#include "matrix3x3.h"
#include "utils.h"
#include "vector3.h"

class PRV3 {
    Vector3 m_axis {};
    double m_principalAngle {};

public:
    static constexpr int SIZE {3};

    PRV3 (const Vector3& axis, const double angleRad) : m_axis {axis}, m_principalAngle {angleRad} {}

    explicit PRV3 (Matrix3x3& mat) {
        //No orthonormality check here, yet. Assumes it is a proper DCM, I do not handle the singularities here.
        m_principalAngle = std::acos((mat.trace()-1)/2);
        m_axis = {
            mat[2, 3] - mat[3, 2],
            mat[3, 1] - mat[1, 3],
            mat[1, 2] - mat[2, 1]
        };

        m_axis = m_axis / (2*std::sin(m_principalAngle));
    }

    friend PRV3 operator+(const PRV3& p1, const PRV3 p2) {
        using std::acos, std::sin, std::cos;

        const double addedAngle {2*acos((cos(p1.m_principalAngle/2)*cos(p2.m_principalAngle/2))
                            -(sin(p1.m_principalAngle/2)*sin(p2.m_principalAngle/2)*(p1.m_axis*p2.m_axis)))};

        if ( std::abs(sin(addedAngle/2)) < 0.01) throw std::invalid_argument("Error: Added angle hits a singularity in PRV.");

        const Vector3 newAxis = (((cos(p2.m_principalAngle/2)*sin(p1.m_principalAngle/2))*p1.m_axis)+
                            ((cos(p1.m_principalAngle/2)*sin(p2.m_principalAngle/2))*p2.m_axis)+
                                ((sin(p2.m_principalAngle/2)*sin(p1.m_principalAngle/2))*(p1.m_axis%p2.m_axis)))/sin(addedAngle/2);

        return PRV3 {newAxis, addedAngle};
    }

    friend std::ostream& operator<<(std::ostream& out, const PRV3& prv3) {
        out << "Axis: {";

        bool hasComma {false};
        for (int i = 0; i < SIZE; ++i) {
            if (hasComma) out << ", ";

            out << prv3.m_axis[i];
            hasComma = true;
        }

        out << "}, Angle: " << prv3.m_principalAngle << std::endl;
        return out;
    }

    [[nodiscard]] double getAngle() const {return m_principalAngle;}
    [[nodiscard]] const Vector3& getAxis() const {return m_axis;}
};
