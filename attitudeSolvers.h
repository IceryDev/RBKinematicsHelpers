#pragma once

#include <assert.h>
#include <iostream>

#include "matrix.h"
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

    inline Matrix getDevenportKMatrix(const std::vector<double>& weights, const std::vector<Vector3>& bodyFrameVectors, const std::vector<Vector3>& inertialVectors) {
        assert(weights.size() == bodyFrameVectors.size() && bodyFrameVectors.size() == inertialVectors.size());

        Matrix3x3 bMatrix {};
        for (std::size_t i {0}; i < weights.size(); ++i) {
            bMatrix += weights.at(i)*(Matrix3x3::outerProduct(bodyFrameVectors.at(i), inertialVectors.at(i)));
        }

        const double sigma {bMatrix.trace()};
        Matrix3x3 sMatrix {(bMatrix+Matrix3x3::transpose(bMatrix))-(sigma*Matrix3x3::identity())};
        Vector3 zMatrix {bMatrix[2, 3]-bMatrix[3, 2], bMatrix[3, 1]-bMatrix[1, 3], bMatrix[1, 2]-bMatrix[2, 1]};

        return Matrix {{
            {sigma, zMatrix[0], zMatrix[1], zMatrix[2]},
            {zMatrix[0], sMatrix[1, 1], sMatrix[1, 2], sMatrix[1, 3]},
            {zMatrix[1], sMatrix[2, 1], sMatrix[2, 2], sMatrix[2, 3]},
            {zMatrix[2], sMatrix[3, 1], sMatrix[3, 2], sMatrix[3, 3]}
        }};
    }
}
