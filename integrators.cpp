#include <vector>
#include <cmath>
#include "integrators.h"
#include "matrix3x3.h"
#include "utils.h"
using namespace Utils;

namespace Integrators {
    std::vector<double> eulerAngles321 {deg2Rad(40), deg2Rad(30), deg2Rad(80)};

    std::vector<double> getTransformedVector(const double time, const double step) {
        using std::sin, std::cos;
        const double radTheta {eulerAngles321.at(1)};
        const double radAlpha {eulerAngles321.at(2)};

        Matrix3x3 m { 0, sin(radAlpha), cos(radAlpha),
                    0, cos(radAlpha) * cos(radTheta), (-sin(radAlpha)) * cos(radTheta),
                    cos(radTheta), sin(radAlpha) * sin(radTheta), cos(radAlpha) * sin(radTheta)};
        std::vector angularVelB {sin(0.1 * time), 0.01, cos(0.1 * time)};

        angularVelB = m.vectorTransform(angularVelB);
        for (auto& angular_vel_b: angularVelB) {
            angular_vel_b *= (deg2Rad(20)/cos(radTheta)) * step;
        }

        return angularVelB;
    }

}