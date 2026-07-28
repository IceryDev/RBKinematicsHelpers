#pragma once
#include "matrix3x3.h"

namespace Integrators {
    extern std::vector<double> eulerAngles321;

    std::vector<double> getTransformedVector(double time, double step);

    inline auto addToState = [] (std::vector<double>& vals, const std::vector<double>& addend) -> std::vector<double>& {
        vals.at(0) += addend.at(0);
        vals.at(1) += addend.at(1);
        vals.at(2) += addend.at(2);

        return vals;
    };

    inline auto targetFunc = [] (const std::vector<double>& vals) {
        return std::__hypot3(vals.at(0), vals.at(1), vals.at(2));
    };

    inline auto integrator(const double duration, const double step) {
        double t {0};
        while (t < duration) {
            addToState(eulerAngles321, getTransformedVector(t, step));
            t += step;
        }

        return targetFunc(eulerAngles321);
    }
}