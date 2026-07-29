#pragma once
#include <cmath>
#include <stdexcept>
#include <initializer_list>
#include <utility>

#include "PRV.h"
#include "utils.h"

class Quaternion {

    class QMatrix {
        static constexpr int DIMENSION {4};

        static constexpr std::array<int, 16> LEFT_MAT {
            0, -1, -2, -3,
            1,  0, -3,  2,
            2,  3,  0, -1,
            3, -2,  1,  0
        };

        static constexpr std::array<int, 16> RIGHT_MAT {
            0, -1, -2, -3,
            1,  0,  3, -2,
            2, -3,  0,  1,
            3,  2, -1,  0
        };

        std::vector<double> m_innerMatrix {};

    public:
        QMatrix(const Quaternion& q, const bool direction) { //true for left, false for right, handled internally
            auto& usedMat {(direction) ? LEFT_MAT : RIGHT_MAT};
            m_innerMatrix.resize(DIMENSION*DIMENSION);

            for (int i = 0; i < DIMENSION*DIMENSION; ++i) {
                m_innerMatrix.at(i) = ((usedMat.at(i) >= 0) ? 1 : -1) * q[std::abs(usedMat.at(i))];
            }
        }

        double& operator[](const int row, const int col) {
            return m_innerMatrix.at((row-1)*DIMENSION+(col-1));
        }
    };

    double m_q0 {};
    double m_q1 {};
    double m_q2 {};
    double m_q3 {};

    bool m_optShortAngle {true};
    bool m_isUnit {true};

public:
    static constexpr int PARAM_COUNT {4};

    Quaternion() = default;

    Quaternion(const double q0, const double q1, const double q2, const double q3) :
        m_q0(q0), m_q1(q1), m_q2(q2), m_q3(q3) {
        if (std::abs(Utils::hypot4(q0, q1, q2, q3)-1) > 0.01) m_isUnit = false;
    }

    Quaternion(const std::initializer_list<double> list) {
        if (list.size() > PARAM_COUNT) throw std::range_error("Error: Quaternion cannot have more than 4 parameters!");

        for (int i = 0; i < list.size(); ++i) {
            (*this)[i] = list.begin()[i];
        }

        if (std::abs(Utils::hypot4(m_q0, m_q1, m_q2, m_q3)-1) > 0.01) m_isUnit = false;
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

        maxVal = sqrt(maxVal);
        (*this)[maxIndex] = maxVal;

        for (int i = 0; i < PARAM_COUNT; ++i) {
            if (i == maxIndex) continue;
            (*this)[i] /= 4*maxVal;
        }

        if ((*this)[0] < 0 && m_optShortAngle) (*this) = -(*this);
        if ((*this)[0] > 0 && !m_optShortAngle) (*this) = -(*this);

        //std::cout << Utils::hypot4(m_q0, m_q1, m_q2, m_q3) << '\n'; //Debug
    }

    void validate() {
        m_isUnit = (std::abs(Utils::hypot4(m_q0, m_q1, m_q2, m_q3)-1) <= 1e-7);
    }

    [[nodiscard]] Matrix3x3 getDCM() const {
        using std::pow;
        if (!m_isUnit) throw std::invalid_argument("Error: A non-unit Quaternion does not describe attitude, therefore DCM cannot be generated!");

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

    const double& operator[](const int index) const {
        switch (index) {
            case 0: return m_q0;
            case 1: return m_q1;
            case 2: return m_q2;
            case 3: return m_q3;
            default: throw std::range_error("Error: Quaternion indexing out of bounds!");
        }
    }

    double& operator[](const int index) {
        m_isUnit = false;
        return const_cast<double&>(std::as_const(*this)[index]);
    }

    Quaternion operator!() const {
        return Quaternion {m_q0, -m_q1, -m_q2, -m_q3};
    }

    Quaternion operator-() const {
        return Quaternion {-m_q0, -m_q1, -m_q2, -m_q3};
    }

    //Keep in mind the rotations are applied from left to right, so q1 first then q2 would be q2*q1
    //Which would have a DCM equivalent of [FN(q)]=[FB(q2)][BN(q1)]
    friend Quaternion operator*(const Quaternion& q1, const Quaternion& q2) {
        Quaternion result {};
        QMatrix leftM {q1, true};
        for (int r = 1; r <= 4; r++){
            double cell = 0;
            for (int i = 1; i <= 4; i++) {
                cell += leftM[r, i] * q2[i-1];
            }
            result[r-1] = cell;
        }

        if (result[0] < 0 && (q1.m_optShortAngle || q2.m_optShortAngle)) result = -result;
        if (result[0] > 0 && (!q1.m_optShortAngle && !q2.m_optShortAngle)) result = -result;

        result.m_isUnit = (!q1.m_isUnit || !q2.m_isUnit);
        return result;
    }

    friend Quaternion operator*(const Quaternion& q, const double val) {
        return Quaternion {q[0]*val, q[1]*val, q[2]*val, q[3]*val};
    }

    friend Quaternion operator/(const Quaternion& q, const double val) {
        if (std::abs(val) < 1e-6) throw std::invalid_argument("Error: Cannot divide by zero!");

        return Quaternion {q[0]/val, q[1]/val, q[2]/val, q[3]/val};
    }

    friend std::ostream& operator<<(std::ostream& out, const Quaternion& q) {
        out << std::format("Quaternion [{}, {}, {}, {}]", q.m_q0, q.m_q1, q.m_q2, q.m_q3);
        return out;
    }

    Quaternion& setShortAngleOption(const bool val) {
        m_optShortAngle = val;
        return (*this);
    }

    [[nodiscard]] QMatrix getLeftMatrix() const {
        return QMatrix {(*this), true};
    }

    //This is only for testing purposes
    [[nodiscard]] Quaternion _applyDKE(const Vector3& w) const {
        return ((*this) * Quaternion {0, w[0], w[1], w[2]})/2;
    }

};
