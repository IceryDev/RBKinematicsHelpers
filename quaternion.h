#pragma once
#include <cmath>
#include <stdexcept>
#include <initializer_list>

#include "PRV.h"
#include "utils.h"

class Quaternion {

    double m_q0 {};
    double m_q1 {};
    double m_q2 {};
    double m_q3 {};

public:
    static constexpr int PARAM_COUNT {4};

    Quaternion(const double q0, const double q1, const double q2, const double q3) :
        m_q0(q0), m_q1(q1), m_q2(q2), m_q3(q3) {
        if (std::abs(Utils::hypot4(q0, q1, q2, q3)-1) > 0.01) throw std::invalid_argument("Error: Quaternion parameters do not satisfy quaternion constraints!");
    }

    Quaternion(const std::initializer_list<double> list) {
        if (list.size() > PARAM_COUNT) throw std::range_error("Error: Quaternion cannot have more than 4 parameters!");

        for (int i = 0; i < list.size(); ++i) {
            (*this)[i] = list.begin()[i];
        }

        if (std::abs(Utils::hypot4(m_q0, m_q1, m_q2, m_q3)-1) > 0.01) throw std::invalid_argument("Error: Quaternion parameters do not satisfy quaternion constraints!");
    }

    explicit Quaternion(const PRV3& prv3) {
        using std::cos, std::sin;
        const double angle { prv3.getAngle() };
        m_q0 = std::cos(angle/2);
        m_q1 = prv3.getAxis()[0]*std::sin(angle/2);
        m_q2 = prv3.getAxis()[1]*std::sin(angle/2);
        m_q3 = prv3.getAxis()[2]*std::sin(angle/2);
    }

    explicit Quaternion(Matrix3x3& mat) {
        //Uses Sheppard's Method, assumes mat is a proper DCM
        using std::sqrt;

        const double trace { mat.trace() };

        //This is step 1: get farthest from zero
        int maxIndex {0};
        double maxVal {(1+trace)/4};
        for (int i = 1; i < PARAM_COUNT; ++i) {
            if (const double tmp {(1+2*mat[i, i]-trace)/4}; maxVal < tmp) {
                maxIndex = i;
                maxVal = tmp;
            }
        }

        auto q0q1 = [=] () {
            return (mat[2, 3]-mat[3, 2]);
        };
        auto q0q2 = [=] () {
            return (mat[3, 1]-mat[1, 3]);
        };
        auto q0q3 = [=] () {
            return (mat[1, 2]-mat[2, 1]);
        };
        auto q1q2 = [=] () {
            return (mat[1, 2]+mat[2, 1]);
        };
        auto q3q1 = [=] () {
            return (mat[3, 1]+mat[1, 3]);
        };
        auto q2q3 = [=] () {
            return (mat[2, 3]+mat[3, 2]);
        };

        maxVal = sqrt(maxVal);
        switch (maxIndex) {
            case 0:
                (*this)[1] = q0q1();
                (*this)[2] = q0q2();
                (*this)[3] = q0q3();
                break;
            case 1:
                (*this)[0] = q0q1();
                (*this)[2] = q1q2();
                (*this)[3] = q3q1();
                break;
            case 2:
                (*this)[0] = q0q2();
                (*this)[1] = q1q2();
                (*this)[3] = q2q3();
                break;
            case 3:
                (*this)[0] = q0q3();
                (*this)[1] = q3q1();
                (*this)[2] = q2q3();
                break;
            default:
                throw std::logic_error("Error: This wasn't supposed to happen! Check your logic!");
        }

        for (int i = 0; i < PARAM_COUNT; ++i) {
            (*this)[i] /= 4*(maxVal);
        }
    }

    [[nodiscard]] Matrix3x3 getDCM() const {
        using std::pow;
        return Matrix3x3 {
            pow(m_q0, 2)+pow(m_q1, 2)-pow(m_q2, 2)-pow(m_q3, 2),
            2*(m_q1*m_q2+m_q0*m_q3),
            2*(m_q1*m_q3-m_q0*m_q2),
            2*(m_q1*m_q2-m_q0*m_q3),
            pow(m_q0, 2)-pow(m_q1, 2)+pow(m_q2, 2)-pow(m_q3, 2),
            2*(m_q2*m_q3+m_q0*m_q1),
            2*(m_q1*m_q3+m_q0*m_q2),
            2*(m_q2*m_q3-m_q0*m_q1),
            pow(m_q0, 2)-pow(m_q1, 2)-pow(m_q2, 2)+pow(m_q3, 2)
        };
    }

    explicit operator Matrix3x3() const {
        return getDCM();
    }

    double& operator[](const int index) {
        switch (index) {
            case 0: return m_q0;
            case 1: return m_q1;
            case 2: return m_q2;
            case 3: return m_q3;
            default: throw std::range_error("Error: Quaternion indexing out of bounds!");
        }
    }

    friend std::ostream& operator<<(std::ostream& out, const Quaternion q) {
        out << std::format("Quaternion [{}, {}, {}, {}]\n", q.m_q0, q.m_q1, q.m_q2, q.m_q3);
        return out;
    }

};
