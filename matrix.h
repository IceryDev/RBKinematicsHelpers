#pragma once
#include <initializer_list>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <format>

class Matrix {
    std::size_t m_rows {};
    std::size_t m_cols {};

    std::vector<double> m_arr {};

public:
    Matrix(const std::size_t row, const std::size_t col) : m_arr(row*col), m_cols(col), m_rows(row) {}

    Matrix(const std::initializer_list<std::initializer_list<double>> elements) {
        m_rows = elements.size();
        m_cols = elements.begin()->size();

        m_arr.resize(m_cols*m_rows);
        auto currentPlace {m_arr.begin()};
        for (const auto element : elements) {
            if (element.size() != m_cols) throw std::invalid_argument("Error: Matrix cannot have rows/columns of different sizes!");
            std::ranges::copy(element, currentPlace);
            currentPlace += static_cast<long>(m_cols);
        }
    }

    Matrix matMul(Matrix matrix){

        if (this->m_cols != matrix.m_rows) throw std::invalid_argument("Error: Invalid matrix dimensions for matrix multiplication!");
        Matrix tempMatrix {m_rows, matrix.m_cols};
        for (int r = 1; r <= m_rows; r++){
            for (int c = 1; c <= matrix.m_cols; c++){
                double cell = 0;
                for (int i = 1; i <= matrix.m_rows; i++) {
                    cell += (*this)[r, i] * matrix[i, c];
                }
                tempMatrix[r, c] = cell;
            }
        }

        return tempMatrix;
    }

    double& operator[](const int row, const int col) {
        return m_arr.at((row-1)*m_cols+(col-1));
    }

    double operator[](const int row, const int col) const {
        return m_arr.at((row-1)*m_cols+(col-1));
    }

    std::pair<std::size_t, std::size_t> dim() {
        return {m_rows, m_cols};
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& m) {

        out << std::format("Matrix{}x{}: [\n", m.m_rows, m.m_cols);
        for (int i {1}; i <= m.m_rows; ++i) {
            out << '\t';
            bool comma {false};
            for (int j {1}; j <= m.m_cols; ++j) {
                if (comma) out << ", ";
                out << m[i, j];
                comma = true;
            }
            out << '\n';
        }
        out << ']';

        return out;
    }
};
