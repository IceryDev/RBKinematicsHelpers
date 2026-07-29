#pragma once
#include <numbers>
#include <vector>
#include <cmath>
#include <stdexcept>

namespace Utils {
    inline auto deg2Rad = [] (const double deg) {return deg * std::numbers::pi/180;};
    inline auto rad2Deg = [] (const double rad) {return (rad / std::numbers::pi) * 180;};

    inline std::vector<double> getCrossProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
        if (v1.size() != 3 || v2.size() != 3) throw std::invalid_argument("Error: Cross product is defined for 3D vectors!");

        return {v1.at(1)*v2.at(2)-v1.at(2)*v2.at(1),
                   v1.at(2)*v2.at(0)-v1.at(0)*v2.at(2),
                   v1.at(0)*v2.at(1)-v1.at(1)*v2.at(0)};
    }

    inline double getDotProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
        if (v1.size() != v2.size()) throw std::invalid_argument("Error: Dot product cannot be done to vectors with different dimensions!");

        double total {0};
        for (std::size_t i = 0; i < v1.size(); ++i) {
            total += v1.at(i)*v2.at(i);
        }

        return total;
    }

    inline double hypot4(const double a1, const double a2, const double a3, const double a4) {
        return std::sqrt(std::pow(a1, 2)+std::pow(a2, 2) + std::pow(a3, 2) + std::pow(a4, 2));
    }
}