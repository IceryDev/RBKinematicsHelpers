#pragma once
#include <initializer_list>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <format>

#include "matrix3x3.h"

class Matrix {

    class MRow {
        Matrix& m_refMatrix;
        std::size_t m_row;

    public:
        MRow(Matrix& mat, const std::size_t row) : m_refMatrix(mat), m_row(row){
            if (mat.m_rows < row) throw std::out_of_range(std::format("Error: Row number {} is out of range for matrix with {} row(s).", row, mat.m_rows));
        }

        MRow& addMultiple(const double scalar, const MRow& mr) {
            if (mr.m_refMatrix.m_cols != m_refMatrix.m_cols) throw std::invalid_argument("Error: Cannot add rows with different sizes!");
            for (int i {1}; i <= m_refMatrix.m_cols; ++i) {
                m_refMatrix[static_cast<int>(m_row), i] += (scalar*mr.m_refMatrix[static_cast<int>(mr.m_row), i]);
            }
            return *this;
        }

        friend Matrix operator*(const double scalar, const MRow& mr) {
            Matrix m {mr};
            for (int i {1}; i <= m.m_cols; ++i) {
                m[1, i] *= scalar;
            }
            return m;
        }

        friend Matrix operator/(const MRow& mr, const double scalar) {
            Matrix m {mr};
            for (int i {1}; i <= m.m_cols; ++i) {
                m[1, i] /= scalar;
            }
            return m;
        }

        MRow& operator*=(const double scalar) {
            if (m_refMatrix.m_det) m_refMatrix.m_det.value() *= scalar;
            for (int i {1}; i <= m_refMatrix.m_cols; ++i) {
                m_refMatrix[static_cast<int>(m_row), i] *= scalar;
            }

            return *this;
        }

        MRow& operator/=(const double scalar) {
            if (m_refMatrix.m_det) m_refMatrix.m_det.value() /= scalar;
            for (int i {1}; i <= m_refMatrix.m_cols; ++i) {
                m_refMatrix[static_cast<int>(m_row), i] /= scalar;
            }

            return *this;
        }

        MRow& operator+=(const MRow& mr) {
            if (mr.m_refMatrix.m_cols != m_refMatrix.m_cols) throw std::invalid_argument("Error: Cannot add rows with different sizes!");

            for (int i {1}; i <= m_refMatrix.m_cols; ++i) {
                m_refMatrix[static_cast<int>(m_row), i] += mr.m_refMatrix[static_cast<int>(mr.m_row), i];
            }

            return *this;
        }

        MRow& operator+=(const Matrix& m) {
            if (m.m_cols != m_refMatrix.m_cols) throw std::invalid_argument("Error: Cannot add rows with different sizes!");

            for (int i {1}; i <= m_refMatrix.m_cols; ++i) {
                m_refMatrix[static_cast<int>(m_row), i] += m[1, i];
            }

            return *this;
        }

        [[nodiscard]] std::size_t getWidth() const {
            return m_refMatrix.m_cols;
        }

        friend class Matrix;
    };

    std::size_t m_rows {};
    std::size_t m_cols {};

    std::vector<double> m_arr {};
    std::optional<double> m_det {std::nullopt};

public:
    explicit Matrix() : m_rows(1), m_cols(1), m_arr(1) {}

    Matrix(const std::size_t row, const std::size_t col) : m_rows(row), m_cols(col), m_arr(row*col) {}

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

    explicit Matrix(const Matrix3x3& mat3) : m_rows(3), m_cols(3), m_arr(9){
        std::ranges::copy(mat3.getInnerArray(), m_arr.begin());
    }

    explicit Matrix(const Vector3& v3, const bool colVec = true) : m_rows((colVec) ? 3 : 1), m_cols((colVec) ? 1 : 3), m_arr(3){
        m_arr[0] = v3[0];
        m_arr[1] = v3[1];
        m_arr[2] = v3[2];
    }

    explicit Matrix(const MRow& mr) : m_rows(1), m_cols(mr.getWidth()), m_arr(m_cols*m_rows){
        std::ranges::copy(&mr.m_refMatrix[static_cast<int>(mr.m_row), 1], &mr.m_refMatrix[static_cast<int>(mr.m_row), static_cast<int>(m_cols)] + 1, m_arr.begin());
    }

    [[nodiscard]] std::size_t getCols() const {return m_cols;}
    [[nodiscard]] std::size_t getRows() const {return m_rows;}

    static Matrix identity(const int size) {
        Matrix m {static_cast<std::size_t>(size), static_cast<std::size_t>(size)};

        for (int i {1}; i <= size; ++i) {
            m[i, i] = 1;
        }

        return m;
    }

    void swapRows(const int r1, const int r2) {
        if (m_det) m_det.value() *= -1;
        for (int col {1}; col <= static_cast<int>(m_cols); ++col) {
            std::swap((*this)[r1, col], (*this)[r2, col]);
        }
    }

    [[nodiscard]] double trace() const {
        if (m_rows != m_cols) throw std::invalid_argument("Error: A matrix must be square to get its trace!");
        double sum {0};
        for (int i {1}; i <= m_rows; ++i) sum += (*this)[i, i];
        return sum;
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

    [[nodiscard]] double det() const {
        if (m_rows != m_cols) throw std::invalid_argument("Error: The determinant is only defined for square matrices!");
        if (m_det) return m_det.value();

        Matrix m {*this};
        return getRREF(m, std::nullopt).det();
    }

    Matrix inv() {
        if (m_rows != m_cols) throw std::invalid_argument("Error: The inverse is only defined for square matrices!");
        if (m_rows == 1) return *this;

        Matrix m {*this};
        Matrix inv {identity(static_cast<int>(m_rows))};
        return getRREF(m, inv);
    }

    //Applies row operations when reducing source to target, source will be reduced and target can be used to get inverse of a matrix for instance
    static Matrix& getRREF(Matrix& source, const std::optional<std::reference_wrapper<Matrix>> target = std::nullopt) {
        Matrix& t = (target) ? target.value().get() : source;
        assert(t.m_rows == source.m_rows);

        reduceRow(static_cast<int>(source.m_rows), source, target);
        backSubstitute(source, target);
        return t;
    }

    double& operator[](const int row, const int col) {
        return m_arr.at((row-1)*m_cols+(col-1));
    }

    double operator[](const int row, const int col) const {
        return m_arr.at((row-1)*m_cols+(col-1));
    }

    MRow operator[](const int row) {
        return MRow {*this, static_cast<std::size_t>(row)};
    }

    std::pair<std::size_t, std::size_t> dim() {
        return {m_rows, m_cols};
    }

    Matrix operator-() {
        Matrix m {m_rows, m_cols};
        for (int i {1}; i <= m_rows; ++i) {
            for (int j {1}; j <= m_cols; ++j) {
                m[i, j] = -((*this)[i, j]);
            }
        }
        return m;
    }

    friend Matrix operator+(const Matrix& m1, const Matrix& m2) {
        assert(m1.m_rows == m2.m_rows && m1.m_cols == m2.m_cols);

        Matrix m {m1.m_rows, m1.m_cols};
        for (int i {1}; i <= m1.m_rows; ++i) {
            for (int j {1}; j <= m1.m_cols; ++j) {
                m[i, j] = m1[i, j] + m2[i, j];
            }
        }
        return m;
    }

    friend Matrix operator-(const Matrix& m1, const Matrix& m2) {
        assert(m1.m_rows == m2.m_rows && m1.m_cols == m2.m_cols);

        Matrix m {m1.m_rows, m1.m_cols};
        for (int i {1}; i <= m1.m_rows; ++i) {
            for (int j {1}; j <= m1.m_cols; ++j) {
                m[i, j] = m1[i, j] - m2[i, j];
            }
        }
        return m;
    }

    friend Matrix operator*(const double scalar, const Matrix& m1) {

        Matrix m {m1.m_rows, m1.m_cols};
        for (int i {1}; i <= m1.m_rows; ++i) {
            for (int j {1}; j <= m1.m_cols; ++j) {
                m[i, j] = m1[i, j] * scalar;
            }
        }
        return m;
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

private:
    static void reduceRow(const int rows, Matrix& source, const std::optional<std::reference_wrapper<Matrix>> target = std::nullopt) {

        Matrix* m {(target) ? &target.value().get() : nullptr};
        double det {1};
        for (int row {1}; row <= rows; ++row) {

            if (row != 1) {
                for (int i {1}; i < row; ++i) {
                    const double coef {-source[row, i]};
                    if (m) (*m)[row].addMultiple(coef, (*m)[i]);
                    source[row].addMultiple(coef, source[i]);
                }
            }

            if (source[row, row] == 0) {
                int swapIdx {-1};
                for (int k {row + 1}; k <= rows; ++k) {
                    if (std::abs(source[k, row]) > 1e-10) { swapIdx = k; break; }
                }
                if (swapIdx == -1) {
                    if (m) throw std::invalid_argument("Error: The inputted matrix is singular!");
                    source.m_det = 0.0;
                    return;
                }
                source.swapRows(row, swapIdx);
                det = -det;
                if (m) m->swapRows(row, swapIdx);
            }

            const double coef2 {source[row, row]};
            det *= coef2;
            if (m) (*m)[row] /= coef2;
            source[row] /= coef2;
        }

        if (source.m_rows == source.m_cols) source.m_det = det;
    }

    static void backSubstitute(Matrix& source, const std::optional<std::reference_wrapper<Matrix>> target = std::nullopt) {

        Matrix* m {(target) ? &target.value().get() : nullptr};

        for (int i {1}; i < source.m_rows; ++i) {

            for (int j {i+1}; j <= source.m_rows; ++j) {
                if (m) (*m)[i].addMultiple(-source[i, j], (*m)[j]);
                source[i].addMultiple(-source[i, j], source[j]);
            }
        }
    }
};
