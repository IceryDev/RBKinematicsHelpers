#pragma once

#include <cmath>
#include <iostream>
#include <format>
#include "matrix3x3.h"
#include "vector3.h"

class Quaternion;

class CRP {
    double m_q1 {};
    double m_q2 {};
    double m_q3 {};

public:
    explicit CRP() = default;

    CRP(const double q1, const double q2, const double q3) :
        m_q1(q1), m_q2(q2), m_q3(q3) {}

    explicit CRP(Matrix3x3& mat) {
        const double coef2 {mat.trace()+1};

        if (std::abs(coef2) < 1e-7)
            throw std::invalid_argument("Error: Given DCM hits a singularity of the Classical Rodrigues Parameters! Consider using a Quaternion instead!");

        m_q1 = (mat[2, 3]-mat[3, 2])/coef2;
        m_q2 = (mat[3, 1]-mat[1, 3])/coef2;
        m_q3 = (mat[1, 2]-mat[2, 1])/coef2;
    }

    [[nodiscard]] Matrix3x3 getDCM() const;

    friend std::ostream& operator<<(std::ostream& out, const CRP& crp) {
        out << std::format("CRP: [{}, {}, {}]", crp.m_q1, crp.m_q2, crp.m_q3);
        return out;
    }

    friend CRP operator+(const CRP& c1, const CRP& c2) {
        Vector3 q1 {c1.m_q1, c1.m_q2, c1.m_q3};
        const Vector3 q2 {c2.m_q1, c2.m_q2, c2.m_q3};

        q1 = (q1+q2-(q1%q2))/(1-(q1*q2));

        return CRP {q1[0], q1[1], q1[2]};
    }
};
