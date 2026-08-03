#pragma once
#include <initializer_list>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <format>

class Matrix {

    class MRow {
        Matrix& m_refMatrix;
        std::size_t m_row;

    public:
        MRow(Matrix& mat, const std::size_t row) : m_refMatrix(mat), m_row(row){
            if (mat.m_rows < row) throw std::out_of_range(std::format("Error: Row number {} is out of range for matrix with {} row(s).", row, mat.m_rows));
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
            for (int i {1}; i <= m_refMatrix.m_cols; ++i) {
                m_refMatrix[static_cast<int>(m_row), i] *= scalar;
            }

            return *this;
        }

        MRow& operator/=(const double scalar) {
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

    explicit Matrix(const MRow& mr) : m_rows(1), m_cols(mr.getWidth()), m_arr(m_cols*m_rows){
        std::ranges::copy(&mr.m_refMatrix[static_cast<int>(mr.m_row), 1], &mr.m_refMatrix[static_cast<int>(mr.m_row), static_cast<int>(m_cols)] + 1, m_arr.begin());
    }

    static Matrix identity(const int size) {
        Matrix m {static_cast<std::size_t>(size), static_cast<std::size_t>(size)};

        for (int i {1}; i <= size; ++i) {
            m[i, i] = 1;
        }

        return m;
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
        for (int row {1}; row <= rows; ++row) {


            if (row != 1) {
                for (int i {1}; i < row; ++i) {
                    if (m) (*m)[row] += -source[row, i] * (*m)[i];
                    source[row] += -source[row, i] * source[i];
                }
            }

            for (int k {1}; k < source.m_rows; ++k) {
                if (source[row, row] != 0) break;
                int idx = ((row - 1 + k) % source.m_rows) + 1;
                source[row] += source[idx];
                if (m) (*m)[row] += (*m)[idx];
            }

            if (source[row, row] == 0) {
                throw std::invalid_argument("Error: The inputted matrix is singular!");
                return;
            }

            if (m) (*m)[row] /= source[row, row];
            source[row] /= source[row, row];
        }
    }

    static void backSubstitute(Matrix& source, const std::optional<std::reference_wrapper<Matrix>> target = std::nullopt) {

        Matrix* m {(target) ? &target.value().get() : nullptr};

        for (int i {1}; i < source.m_rows; ++i) {
            for (int j {i+1}; j <= source.m_rows; ++j) {
                if (m) (*m)[i] += -source[i, j] * (*m)[j];
                source[i] += -source[i, j] * source[j];
            }
        }
    }
};
