#include <string_view>
#include <sstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>


#include "lcodes.h"

namespace coding {

    namespace gsl_wrapper {
        gsl_matrix_ptr gsl_matrix_transpose_non_square(const gsl_matrix &mx) {
            coding::gsl_matrix_ptr transposed{
                    gsl_matrix_alloc(mx.size2, mx.size1),
                    &gsl_matrix_free
            };

            for (size_t i = 0; i < mx.size2; i++)
                for (size_t j = 0; j < mx.size1; j++)
                    gsl_matrix_set(transposed.get(), i, j, gsl_matrix_get(&mx, j, i));

            return transposed;
        }
    }

    namespace util {

        // Requires generator_mx to be in standard form: G = (E|A).
        // Given a generator_mx matrix with dimensions kxn, the resulting parity_mx check matrix is tx(k+t) where t := n - k.

        gsl_matrix_ptr generator_to_parity(const gsl_matrix &generator_mx) {
            const std::size_t k = generator_mx.size1;
            const std::size_t n = generator_mx.size2;
            const std::size_t t = n - k;

            coding::gsl_matrix_ptr parity_mx{
                    gsl_matrix_alloc(t, k + t),
                    &gsl_matrix_free
            };

            gsl_matrix_const_view generator_A_cview = gsl_matrix_const_submatrix(&generator_mx, 0, k, k, n - k);
            gsl_matrix_ptr temp_A{gsl_matrix_alloc(k, t), &gsl_matrix_free};
            gsl_matrix_memcpy(temp_A.get(), &generator_A_cview.matrix);
            auto At = gsl_wrapper::gsl_matrix_transpose_non_square(*temp_A);

            gsl_matrix_view parity_A_view = gsl_matrix_submatrix(parity_mx.get(), 0, 0, t, k);
            gsl_matrix_memcpy(&parity_A_view.matrix, At.get());

            gsl_matrix_view parity_E_view = gsl_matrix_submatrix(parity_mx.get(), 0, k, t, t);
            gsl_matrix_set_identity(&parity_E_view.matrix);

            return parity_mx;
        }

        gsl_matrix_ptr parity_to_generator(const gsl_matrix &parity_mx) {
            const std::size_t t = parity_mx.size1;
            const std::size_t k = parity_mx.size2 - parity_mx.size1;
            const std::size_t n = t + k;

            coding::gsl_matrix_ptr generator_mx{
                    gsl_matrix_alloc(k, n),
                    &gsl_matrix_free
            };

            gsl_matrix_const_view parity_A_cview = gsl_matrix_const_submatrix(&parity_mx, 0, 0, t, k);
            gsl_matrix_ptr temp_A{gsl_matrix_alloc(t, k), &gsl_matrix_free};
            gsl_matrix_memcpy(temp_A.get(), &parity_A_cview.matrix);
            auto At = gsl_wrapper::gsl_matrix_transpose_non_square(*temp_A);

            gsl_matrix_view generator_A_view = gsl_matrix_submatrix(generator_mx.get(), 0, k, k, t);
            gsl_matrix_memcpy(&generator_A_view.matrix, At.get());

            gsl_matrix_view generator_E_view = gsl_matrix_submatrix(generator_mx.get(), 0, 0, k , k);
            gsl_matrix_set_identity(&generator_E_view.matrix);

            return generator_mx;
        }

        gsl_matrix_ptr generator_to_standard_form(gsl_matrix_ptr generator) {
            return generator;
        }
    }

    namespace gsl_wrapper {

        gsl_vector_ptr gsl_matrix_mul_vector(const gsl_vector &vec, const gsl_matrix &matrix) {
            constexpr static double alpha = 1.0; // scaling factor for the matrix-vector product
            constexpr static double beta = 0.0; // scaling factor for the result vector
            gsl_vector_ptr result{gsl_vector_alloc(matrix.size2), &gsl_vector_free};
            gsl_blas_dgemv(CblasNoTrans, alpha, &matrix, &vec, beta, result.get());
            return result;
        }
    }

    /*
     * Helpers
     */
    void linear_code::evaluate_min_distance() {}

    void linear_code::evaluate_capabilities() {}

    void linear_code::evaluate_decoding_table() {}

    void linear_code::sanity_check_properties() const {
        const bool sane_generator_matrix_dimensions = m_code->size1 > m_code->size2 && m_code->size1 < (1 << 7);
        if (!sane_generator_matrix_dimensions)
            throw linear_code_exception{fmt::format("Error: Insane generator matrix dimensions!")};
    }

    gsl_vector_ptr linear_code::with_redundancy(const gsl_vector &word) const {
        if (word.size != m_code->size1)
            throw linear_code_exception{fmt::format(
                    "Error: Cannot add redundancy to word of size {} using generator matrix {}x{}",
                    word.size, m_code->size1, m_code->size2)};

        auto word_with_redundancy = gsl_wrapper::gsl_matrix_mul_vector(word, *m_code);
        return word_with_redundancy;
    }

    /*
     * Instance control
     */

    class linear_code linear_code::from_parity_equations(gsl_matrix_ptr &&parity) {
        auto generator = util::parity_to_generator(*parity);
        auto code = linear_code{std::move(generator)};
        // Initialize properties.
        return code;
    }

    class linear_code linear_code::from_generator(gsl_matrix_ptr &&generator) {
        auto code = linear_code{std::move(generator)};
        // Initialize properties.
        return code;
    }

    class linear_code linear_code::from_dual(const linear_code &dual) {
    }

    linear_code::linear_code(gsl_matrix_ptr &&generator)
            : m_code{util::generator_to_standard_form(std::move(generator))}, m_basis_size{m_code->size1},
              m_word_length{m_code->size2} {
        // NB: m_code is initialized first because it is declared first in the header file and not
        // 	   because it is put first in the list here.
    }

    /*
     * Operations
     */

    gsl_vector_ptr linear_code::encode(const gsl_vector &word) {
        return with_redundancy(word);
    }

    gsl_vector_ptr linear_code::decode(const gsl_vector &) {}
}
