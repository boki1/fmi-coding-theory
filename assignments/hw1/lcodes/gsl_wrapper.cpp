#include "gsl_wrapper.h"

namespace gsl_wrapper {

    gsl_matrix_ptr gsl_matrix_transpose_non_square(const gsl_matrix &mx) {
        gsl_matrix_ptr transposed{
                gsl_matrix_alloc(mx.size2, mx.size1),
                &gsl_matrix_free
        };

        for (size_t i = 0; i < mx.size2; i++)
            for (size_t j = 0; j < mx.size1; j++)
                gsl_matrix_set(transposed.get(), i, j, gsl_matrix_get(&mx, j, i));

        return transposed;
    }

    gsl_vector_ptr gsl_matrix_mul_vector(const gsl_vector &vec, const gsl_matrix &matrix) {
        assert(vec.size == matrix.size1);
        gsl_vector_ptr result{gsl_vector_alloc(matrix.size2), &gsl_vector_free};

        // This great library does not have an additional matrix-vector multiplication function. The only way I
        // found was to *allocated* a matrix anew and *copy* the vector content there. Then the result will be yet
        // another matrix - whereas it is actually a vector-column, which I am forced to *copy again* into the vector
        // that I want to return. No, thanks.
        for (int j = 0; j < matrix.size2; j++) {
            double sum = 0;
            for (int i = 0; i < matrix.size1; i++)
                sum += gsl_matrix_get(&matrix, i, j) * gsl_vector_get(&vec, i);
            gsl_vector_set(result.get(), j, sum);
        }

        return result;
    }

    // In-place dynamic variant of the other overload.
    void make_word(std::uint64_t bits, gsl_vector &vec) {
        assert(vec.size <= sizeof(bits) * CHAR_BIT);
        for (int i = 0; i < vec.size; ++i)
            gsl_vector_set(&vec, i, bits & (1 << i));
    }


}