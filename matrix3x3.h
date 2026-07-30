#pragma once
#include <vector>
#include <array>
#include <iostream>

#include "vector3.h"

class Matrix3x3 {
    static constexpr int DIMENSION {3};
    std::array<double, DIMENSION*DIMENSION> m_arr {};


public:
    Matrix3x3() = default;

    Matrix3x3(std::initializer_list<double> entries) {
        if (entries.size() > DIMENSION*DIMENSION) throw std::length_error("Error: Matrix3x3 initialised with more than 9 arguments!");

        std::ranges::copy(entries, m_arr.begin());
    }

    static Matrix3x3 identity() {
        return Matrix3x3 {
            1, 0, 0,
            0, 1, 0,
            0, 0, 1
        };
    }

    static Matrix3x3 skewSymmetric(const Vector3& v) {
        return Matrix3x3 {
            0, -v[2], v[1],
            v[2], 0, -v[0],
            -v[1], v[0], 0
        };
    }

    //Calculates [v][transpose(v)], did not use a loop for copy elision
    static Matrix3x3 outerProduct(const Vector3& v) {
        return Matrix3x3 {
            v[0]*v[0], v[0]*v[1], v[0]*v[2],
            v[1]*v[0], v[1]*v[1], v[1]*v[2],
            v[2]*v[0], v[2]*v[1], v[2]*v[2]
        };
    }

    Matrix3x3& transpose() {
        std::swap((*this)[1, 2], (*this)[2, 1]);
        std::swap((*this)[1, 3], (*this)[3, 1]);
        std::swap((*this)[2, 3], (*this)[3, 2]);
        return *this;
    }

    double trace() {
        return (*this)[1, 1] + (*this)[2, 2] + (*this)[3, 3];
    }

    Matrix3x3 matMul(Matrix3x3 matrix){

        Matrix3x3 tempMatrix {};
        for (int r = 1; r <= DIMENSION; r++){
            for (int c = 1; c <= DIMENSION; c++){
                double cell = 0;
                for (int i = 1; i <= DIMENSION; i++) {
                    cell += (*this)[r, i] * matrix[i, c];
                }
                tempMatrix[r, c] = cell;
            }
        }

        return tempMatrix;
    }

    //Keeping this for compatibility with previous Euler Angle integration code
    [[deprecated]] std::vector<double> vectorTransform(const std::vector<double>& vec) {
        std::vector<double> result(DIMENSION);
        for (int r = 1; r <= DIMENSION; r++){
            double cell = 0;
            for (int i = 1; i <= DIMENSION; i++) {
                cell += (*this)[r, i] * vec.at(i-1);
            }
            result.at(r-1) = cell;
        }

        return result;
    }

    Vector3 vectorTransform(const Vector3& vec) {
        Vector3 result {};
        for (int r = 1; r <= DIMENSION; r++){
            double cell = 0;
            for (int i = 1; i <= DIMENSION; i++) {
                cell += (*this)[r, i] * vec[i-1];
            }
            result[r-1] = cell;
        }

        return result;
    }

    double& operator[](const int row, const int col) {
        return m_arr.at((row-1)*DIMENSION+(col-1));
    }

    double operator[](const int row, const int col) const {
        return m_arr.at((row-1)*DIMENSION+(col-1));
    }

    friend Matrix3x3 operator+(const Matrix3x3& m1, const Matrix3x3& m2) {
        //Didn't use a loop for copy elision. Named variables disable elision
        return Matrix3x3 {
            m1[1, 1]+m2[1, 1], m1[1, 2]+m2[1, 2], m1[1, 3]+m2[1, 3],
            m1[2, 1]+m2[2, 1], m1[2, 2]+m2[2, 2], m1[2, 3]+m2[2, 3],
            m1[3, 1]+m2[3, 1], m1[3, 2]+m2[3, 2], m1[3, 3]+m2[3, 3]
        };
    }

    Matrix3x3& divBy(const double scalar) {
        if (std::abs(scalar) < 1e-7) throw std::invalid_argument("Error: Cannot divide by zero!");
        for (auto& val: m_arr) {
            val /= scalar;
        }
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix3x3 mat) {
        for (int i = 0; i < DIMENSION*DIMENSION; ++i) {

            out << mat.m_arr.at(i) << ' ';
            if (i%3 == 2) out << '\n';
        }

        return out;
    }
};
