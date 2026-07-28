#pragma once
#include <vector>
#include <array>
#include <iostream>

class Matrix3x3 {
    static constexpr int DIMENSION {3};
    std::array<double, DIMENSION*DIMENSION> m_arr {};


public:
    Matrix3x3() = default;

    Matrix3x3(std::initializer_list<double> entries) {
        if (entries.size() > DIMENSION*DIMENSION) throw std::length_error("Error: Matrix3x3 initialised with more than 9 arguments!");

        std::ranges::copy(entries, m_arr.begin());
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

    std::vector<double> vectorTransform(const std::vector<double>& vec) {
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

    double& operator[](const int row, const int col) {
        return m_arr.at((row-1)*DIMENSION+(col-1));
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix3x3 mat) {
        for (int i = 0; i < DIMENSION*DIMENSION; ++i) {

            out << mat.m_arr.at(i) << ' ';
            if (i%3 == 2) out << '\n';
        }

        return out;
    }
};