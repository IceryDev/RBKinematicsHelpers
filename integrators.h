#pragma once
#include "matrix3x3.h"
#include "quaternion.h"

namespace Integrators {
    extern std::vector<double> eulerAngles321;
    extern Quaternion eulerParameters;

    std::vector<double> getTransformedVector(double time, double step);

    inline auto addToState = [] (std::vector<double>& vals, const std::vector<double>& addend) -> std::vector<double>& {
        vals.at(0) += addend.at(0);
        vals.at(1) += addend.at(1);
        vals.at(2) += addend.at(2);

        return vals;
    };

    //If you get errors, __hypot3 is implementation defined, that is probably the reason
    inline auto targetFuncEA = [] (const std::vector<double>& vals) {
        return std::__hypot3(vals.at(0), vals.at(1), vals.at(2));
    };

    inline auto targetFunc = [] (const Quaternion& vals) {
        return std::__hypot3(vals[1], vals[2], vals[3]);
    };

    inline auto integrator(const double duration, const double step) {
        double t {0};
        while (t < duration) {
            addToState(eulerParameters, getTransformedVector(t, step));
            t += step;
        }

        return targetFuncEA(eulerAngles321);
    }
}
