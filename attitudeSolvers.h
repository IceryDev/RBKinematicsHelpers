#pragma once

#include <iostream>
#include "matrix3x3.h"
#include "vector3.h"

namespace AttitudeSolvers {
    inline Matrix3x3 triadSolver(Vector3 v1N, Vector3 v1B, Vector3 v2N, Vector3 v2B) {

        v1B.normalise();
        v2B.normalise();
        v1N.normalise();
        v2N.normalise();

        if (std::abs((v1B%v2B).norm()) < 1e-7) throw std::invalid_argument("Error: Two non-parallel observation vectors must be used for a static TRIAD attitude determination.");

        const Vector3 t1B {v1B};
        const Vector3 t2B {(v1B%v2B).normalise()};

        Matrix3x3 matBN {t1B, t2B, (t1B%t2B).normalise()};

        //Inertial frame triad vector
        const Vector3 t1N {v1N};
        const Vector3 t2N {(v1N%v2N).normalise()};

        matBN = matBN.matMul((Matrix3x3 {t1N, t2N, (t1N%t2N).normalise()}).transpose());
        return matBN;
    }


}
