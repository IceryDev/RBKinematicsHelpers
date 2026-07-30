#include <iostream>
#include <functional>
#include <cmath>
#include "matrix3x3.h"
#include "integrators.h"
#include "PRV.h"
#include "quaternion.h"
#include "utils.h"
#include "CRP.h"
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

    CRP c {0.1, 0.2, 0.3};

    std::cout << c.getDCM() << std::endl;

    Matrix3x3 m {0.333333, -0.666667, 0.666667,
                 0.871795, 0.487179, 0.0512821,
                 -0.358974, 0.564103, 0.74359};

    CRP c1 {m};
    std::cout << c1 << std::endl;

    CRP c11 {-0.3,0.3,0.1};
    CRP c12 {-0.1,-0.2,-0.3};

    std::cout << c11 + c12 << std::endl;

}