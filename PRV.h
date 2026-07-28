#pragma once
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>

#include "matrix3x3.h"

class PRV3 {
    std::vector<double> m_axis {};
    double m_principalAngle {};

public:
    PRV3 (std::vector<double> axis, const double angleRad) : m_axis {std::move(axis)}, m_principalAngle {angleRad} {
        if (axis.size() != 3) throw std::range_error("Error: Principal rotation axis must be 3D!");
    }

    explicit PRV3 (Matrix3x3& mat) {
        //No orthonormality check here, yet. Assumes it is a proper DCM
        m_principalAngle = std::acos((mat.trace()-1)/2);
        m_axis = {
            mat[2, 3] - mat[3, 2],
            mat[3, 1] - mat[1, 3],
            mat[1, 2] - mat[2, 1]
        };
        for (auto& val: m_axis) {
            val /= (2*std::sin(m_principalAngle));
        }
    }

    friend std::ostream& operator<<(std::ostream& out, const PRV3& prv3) {
        out << "Axis: {";

        bool hasComma {false};
        for (const auto val: prv3.m_axis) {
            if (hasComma) out << ", ";

            out << val;
            hasComma = true;
        }

        out << "}, Angle: " << prv3.m_principalAngle << std::endl;
        return out;
    }
};
