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
#include "matrix.h"
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

    std::vector<double> weights {1, 1};
    std::vector<Vector3> bodyVec {
        {0.8273, 0.5541, -0.0920},
        {-0.8285, 0.5522, -0.0955}
    };
    std::vector<Vector3> inertialVec {
        {-0.1517, -0.9669, 0.2050},
        {-0.8393, 0.4494, -0.3044}
    };

    Quaternion v {-0.131949323884729, 4.202477253067334, -2.508412840978479, 1};
    v.normalise();
    v.validate();
    std::cout << Quaternion::getDCM(v) << std::endl;
}