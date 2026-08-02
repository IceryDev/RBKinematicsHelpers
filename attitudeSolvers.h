#pragma once

#include <iostream>
#include "matrix3x3.h"
#include "vector3.h"

namespace AttitudeSolvers {
    class TRIADSolver {
        //Inertial frame unit vectors
        Vector3 m_v1N {};
        Vector3 m_v2N {};
        //Body frame unit vectors
        Vector3 m_v1B {};
        Vector3 m_v2B {};

        Matrix3x3 m_matBN {};

    public:

        //Assumes first sensor is more accurate than the second sensor.
        TRIADSolver(const Vector3& v1N, const Vector3& v2N, const Vector3& v1B, const Vector3& v2B) :
            m_v1N(v1N), m_v2N(v2N), m_v1B(v1B), m_v2B(v2B) {

            m_v1B.normalise();
            m_v2B.normalise();
            m_v1N.normalise();
            m_v2N.normalise();

            if (std::abs((m_v1B%m_v2B).norm()) < 1e-7) throw std::invalid_argument("Error: Two non-parallel observation vectors must be used for a static TRIAD attitude determination.");

            const Vector3 t1B {m_v1B};
            const Vector3 t2B {(m_v1B%m_v2B).normalise()};

            m_matBN = Matrix3x3 {t1B, t2B, (t1B%t2B).normalise()};

            //Inertial frame triad vector
            const Vector3 t1N {m_v1N};
            const Vector3 t2N {(m_v1N%m_v2N).normalise()};

            m_matBN = m_matBN.matMul((Matrix3x3 {t1N, t2N, (t1N%t2N).normalise()}).transpose());
        }

        friend std::ostream& operator<<(std::ostream& out, const TRIADSolver& ts) {
            out << ts.m_matBN;
            return out;
        }
    };
}
