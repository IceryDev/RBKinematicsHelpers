#pragma once
#include <numbers>

namespace Utils {
    inline auto deg2Rad = [] (const double deg) {return deg * std::numbers::pi/180;};
    inline auto rad2Deg = [] (const double rad) {return (rad / std::numbers::pi) * 180;};
}