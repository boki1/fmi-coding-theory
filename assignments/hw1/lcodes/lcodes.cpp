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
    }

    /*
     * Helpers
     */
    void linear_code::evaluate_min_distance() {
        const std::size_t num_keywords = 1 << m_basis_size;

        std::size_t min_weight = std::numeric_limits<int>::max();
        gsl_vector_ptr min_weighted_word{gsl_vector_alloc(m_code->size1), &gsl_vector_free};
        gsl_vector_ptr current{gsl_vector_alloc(m_code->size1), &gsl_vector_free};

        // Skip the empty word.
        for (uint64_t i = 1; i < num_keywords; ++i) {
            gsl_wrapper::make_word(i, *current);
            auto word_with_redundancy = with_redundancy(*current);
            if (const auto wt = weight(*word_with_redundancy); wt < min_weight) {
                gsl_vector_memcpy(min_weighted_word.get(), current.get());
                min_weight = wt;
            }
        }

        m_min_distance = min_weight;
    }

    void linear_code::evaluate_capabilities() {
        assert(m_min_distance > 0);
        m_max_errors_detect = m_min_distance - 1;
        m_max_errors_correct = (m_min_distance - 1) / 2;
        m_radius = m_word_length - m_basis_size;
    }

    void linear_code::evaluate_decoding_table() {}

    void linear_code::sanity_check_properties() const {
        const bool sane_generator_matrix_dimensions = m_code->size1 > m_code->size2 && m_code->size1 < (1 << 10);
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
        return code;
    }

    class linear_code linear_code::from_generator(gsl_matrix_ptr &&generator) {
        auto code = linear_code{std::move(generator)};
        return code;
    }

    class linear_code linear_code::from_dual(const linear_code &dual) {
    }

    linear_code::linear_code(gsl_matrix_ptr &&generator)
            : m_code{util::generator_to_standard_form(std::move(generator))}
            , m_basis_size{m_code->size1}
            , m_word_length{m_code->size2} {
        // NB: m_code is initialized first because it is declared first in the header file and not
        // 	   because it is put first in the list here.
        evaluate_min_distance();
        evaluate_capabilities();
    }

    /*
     * Operations
     */

    std::size_t linear_code::weight(const gsl_vector &word) {
        std::size_t weight = 0;
        for (int i = 0; i < word.size; ++i)
            weight += static_cast<std::size_t>(gsl_vector_get(&word, i) != 0);
        return weight;
    }


    gsl_vector_ptr linear_code::encode(const gsl_vector &word) {
        return with_redundancy(word);
    }

    gsl_vector_ptr linear_code::decode(const gsl_vector &) {}
}
