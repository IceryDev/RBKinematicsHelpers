#pragma once
#include "matrix3x3.h"
#include "quaternion.h"

namespace Integrators {
    extern std::vector<double> eulerAngles321;
    extern Quaternion eulerParameters;
    extern Vector3 classicalRodrigues;

    std::vector<double> getTransformedVector(double time, double step);
    Quaternion calculateNextStepQ(double time, double step);
    Vector3 calculateNextStepCRP(double time, double step);

    inline auto addToState = [] (std::vector<double>& vals, const std::vector<double>& addend) -> std::vector<double>& {
        vals.at(0) += addend.at(0);
        vals.at(1) += addend.at(1);
        vals.at(2) += addend.at(2);

        return vals;
    };

    //These three functions below are question specific
    inline auto targetFuncEA = [] (const std::vector<double>& vals) {
        return std::hypot(vals.at(0), vals.at(1), vals.at(2));
    };

    inline auto targetFuncQ = [] (const Quaternion& vals) {
        return std::hypot(vals[1], vals[2], vals[3]);
    };

    inline auto targetFuncCRP = [] (const Vector3& vals) {
        return std::hypot(vals[1], vals[2], vals[0]);
    };

    inline auto integratorCRP(const double duration, const double step) {
        double t {0};
        while (t < duration) {
            classicalRodrigues += calculateNextStepCRP(t, step);
            t += step;
        }

        return targetFuncCRP(classicalRodrigues);
    }

    inline auto integratorQ(const double duration, const double step) {
        double t {0};
        while (t < duration) {
            //Euler Estimation
            eulerParameters = eulerParameters + calculateNextStepQ(t, step);
            //Get rid of errors
            eulerParameters.normalise();
            t += step;
        }

        return targetFuncQ(eulerParameters);
    }

    inline auto integratorEA(const double duration, const double step) {
        double t {0};
        while (t < duration) {
            addToState(eulerAngles321, getTransformedVector(t, step));
            t += step;
        }

        return targetFuncEA(eulerAngles321);
    }
}
