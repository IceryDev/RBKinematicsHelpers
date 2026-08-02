#include <iostream>
#include <functional>
#include <cmath>

#include "attitudeSolvers.h"
#include "matrix3x3.h"
#include "integrators.h"
#include "PRV.h"
#include "quaternion.h"
#include "utils.h"
#include "CRP.h"
#include "MRP.h"
using namespace Utils;


Matrix3x3 get321EulerAngles(const double rad1, const double rad2, const double rad3) {
    using std::sin, std::cos;
    Matrix3x3 mat {};

    mat[1, 1] = cos(rad2) * cos(rad1);
    mat[1, 2] = cos(rad2) * sin(rad1);
    mat[1, 3] = -1 * sin(rad2);
    mat[2, 1] = (sin(rad3) * sin(rad2) * cos(rad1)) - (cos(rad3) * sin(rad1));
    mat[2, 2] = (sin(rad3) * sin(rad2) * sin(rad1)) + (cos(rad3) * cos(rad1));
    mat[2, 3] = sin(rad3) * cos(rad2);
    mat[3, 1] = (cos(rad3) * sin(rad2) * cos(rad1)) + (sin(rad3) * sin(rad1));
    mat[3, 2] = (cos(rad3) * sin(rad2) * sin(rad1)) - (sin(rad3) * cos(rad1));
    mat[3, 3] = cos(rad3) * cos(rad2);
    return mat;
}

Matrix3x3 get321EulerAnglesDeg(const double deg1, const double deg2, const double deg3) {return get321EulerAngles(deg2Rad(deg1), deg2Rad(deg2), deg2Rad(deg3));}


int main() {

    AttitudeSolvers::TRIADSolver ts {
        {-0.1517, -0.9669, 0.2050},
        {-0.8393, 0.4494, -0.3044},
        {0.8273, 0.5541, -0.0920},
        {-0.8285, 0.5522, -0.0955}
    };

    std::cout << ts << std::endl;

    Matrix3x3 estimate {
        0.969846, 0.17101, 0.173648,
        -0.200706, 0.96461, 0.17101,
        -0.138258, -0.200706, 0.969846
    };

    Matrix3x3 trueMat {
        0.963592, 0.187303, 0.190809,
        -0.223042, 0.956645, 0.187303,
        -0.147454, -0.223042, 0.963593
    };

    PRV3 rot {estimate.matMul(Matrix3x3::transpose(trueMat))};
    std::cout << rad2Deg(rot.getAngle());
}