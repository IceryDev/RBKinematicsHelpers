#pragma once
#include <initializer_list>
#include <stdexcept>
#include <vector>

class Vector3 {
    static constexpr int SIZE {3};

    double m_x {};
    double m_y {};
    double m_z {};

public:
    Vector3() = default;

    Vector3(const double x, const double y, const double z) : m_x(x), m_y(y), m_z(z) {}

    Vector3(const std::initializer_list<double> entries) {
        if (entries.size() > 3) throw std::length_error("Error: Vector3 initialised with more than 3 arguments!");

        for (int i = 0; i < entries.size(); ++i) {
            (*this)[i] = entries.begin()[i];
        }
    }

    //Cross product
    friend Vector3 operator%(const Vector3& v1, const Vector3& v2) {

        return {v1[1]*v2[2]-v1[2]*v2[1],
                   v1[2]*v2[0]-v1[0]*v2[2],
                   v1[0]*v2[1]-v1[1]*v2[0]};
    }

    friend double operator*(const Vector3& v1, const Vector3& v2) {

        double total {0};
        for (int i = 0; i < SIZE; ++i) {
            total += v1[i]*v2[i];
        }

        return total;
    }

    friend Vector3 operator+(const Vector3& v1, const Vector3& v2) {
        return Vector3{v1.m_x+v2.m_x, v1.m_y+v2.m_y, v1.m_z+v2.m_z};
    }

    friend Vector3 operator-(const Vector3& v1, const Vector3& v2) {
        return Vector3{v1.m_x-v2.m_x, v1.m_y-v2.m_y, v1.m_z-v2.m_z};
    }

    friend Vector3 operator*(const double scalar, const Vector3& v) {
        return Vector3{v[0]*scalar, v[1]*scalar, v[2]*scalar};
    }

    friend Vector3 operator*(const Vector3& v, const double scalar) {
        return scalar * v;
    }

    friend Vector3 operator/(const Vector3& v, const double scalar) {
        return Vector3{v[0]/scalar, v[1]/scalar, v[2]/scalar};
    }

    double& operator[](const int i) {
        switch (i) {
            case 0: return m_x;
            case 1: return m_y;
            case 2: return m_z;
            default: throw std::range_error("Error: Vector3 index out of range.");
        }
    }

    double operator[](const int i) const {
        switch (i) {
            case 0: return m_x;
            case 1: return m_y;
            case 2: return m_z;
            default: throw std::range_error("Error: Vector3 index out of range.");
        }
    }

    explicit operator std::vector<double>() const {
        return {m_x, m_y, m_z};
    }
};
