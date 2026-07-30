#include "CRP.h"

#include "quaternion.h"

[[nodiscard]] Matrix3x3 CRP::getDCM() const {
    using std::pow;
    const double coef {1+pow(m_q1, 2)+pow(m_q2, 2)+pow(m_q3, 2)};

    Quaternion q {1, m_q1, m_q2, m_q3};
    q.m_isUnit = true; //Override checking behaviour to get DCM

    return Quaternion::getDCM(q).divBy(coef);
}
