#pragma once

#include <cmath>
#include <iostream>
#include <format>
#include "matrix3x3.h"

//Modified Rodrigues Parameters
class MRP {

    double m_s1 {};
    double m_s2 {};
    double m_s3 {};

public:
    explicit MRP() = default;

    MRP (const double s1, const double s2, const double s3) : m_s1(s1), m_s2(s2), m_s3(s3) {}
    explicit MRP (const Vector3& v) : m_s1(v[0]), m_s2(v[1]), m_s3(v[2]) {}

    explicit MRP (Matrix3x3& mat) {
        //Assumes a proper DCM, check orthonormality and det(C)=1 for right handed-ness before.

        double coef {std::sqrt(mat.trace()+1)};
        coef = coef * (coef + 2);

        //Use the Shadow Set if hitting a singularity at 2 pi
        if (std::abs(coef) < 1e-7) {
            m_s1 = m_s2 = m_s3 = 0;
            return;
        }

        m_s1 = (mat[2, 3]-mat[3, 2])/coef;
        m_s2 = (mat[3, 1]-mat[1, 3])/coef;
        m_s3 = (mat[1, 2]-mat[2, 1])/coef;

        //No auto switching to short angle yet
    }

    Matrix3x3 getDCM() {
        Matrix3x3 skewSym {Matrix3x3::skewSymmetric({m_s1, m_s2, m_s3})};
        const Matrix3x3 id3x3 {Matrix3x3::identity()};

        return (id3x3) + ((8*(skewSym.matMul(skewSym)))-((4*(1-norm2()))*skewSym)).divBy(std::pow(1+norm2(), 2));
    }

    [[nodiscard]] double norm2() const {
        using std::pow;
        return pow(m_s1, 2)+pow(m_s2, 2)+pow(m_s3, 2);
    }

    [[nodiscard]] int getZone() const {
        if (const double n {norm2()}; std::abs(n-1) < 1e-7) return 0;
        else if (n > 1) return -1;

        return 1;
    }

    void swapSet() {
        const double normSq {this->norm2()};

        for (int i = 1; i <= 3; ++i) {
            (*this)[i] = (-(*this)[i])/normSq;
        }
    }

    static MRP swapSet(const MRP& m) {
        MRP tmp {m};
        tmp.swapSet();
        return tmp;
    }

    const double& operator[](const int index) const {
        switch (index) {
            case 1: return m_s1;
            case 2: return m_s2;
            case 3: return m_s3;
            default: throw std::range_error("Error: MRP indexing out of bounds!");
        }
    }

    double& operator[](const int index) {
        return const_cast<double&>(std::as_const(*this)[index]);
    }

    MRP operator-() const {
        return MRP {-m_s1, -m_s2, -m_s3};
    }

    friend MRP operator+(const MRP& m1, const MRP& m2) {
        const double norm_m1 {m1.norm2()}; //These are norm squared
        const double norm_m2 {m2.norm2()};

        //Get both as short angles
        const MRP& shortM1 = (m1.getZone() == -1) ? swapSet(m1) : m1;
        const MRP& shortM2 = (m2.getZone() == -1) ? swapSet(m2) : m2;

        //Turn them to vectors to do vector math
        const Vector3 m1Vec {shortM1.m_s1, shortM1.m_s2, shortM1.m_s3};
        const Vector3 m2Vec {shortM2.m_s1, shortM2.m_s2, shortM2.m_s3};

        const double denominator {1+(norm_m1*norm_m2)-(2*(m2Vec*m1Vec))};

        //With two short angles there is still a possibility where we hit a singularity, that is when m1 and m2 both has a principal angle of 180
        //Thus we return the shadow set of 360, which is 0 degrees.
        if (std::abs(denominator) < 1e-7) {
            return MRP {0, 0, 0};
        }

        const Vector3 finalVec = (((1-norm_m2)*m1Vec)+((1-norm_m1)*m2Vec)-(2*(m1Vec%m2Vec)))/denominator;
        MRP result {finalVec};
        if (result.getZone() == -1) result.swapSet();

        return result;
    }

    friend MRP operator-(const MRP& m1, const MRP& m2) {
        return m1 + (-m2);
    }

    friend std::ostream& operator<<(std::ostream& out, const MRP& m) {
        out << std::format("MRP: [{}, {}, {}]", m.m_s1, m.m_s2, m.m_s3);
        return out;
    }

};