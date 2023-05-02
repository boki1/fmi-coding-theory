#include <string_view>
#include <algorithm>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>

#include <fmt/core.h>

#include "lcodes.h"

namespace coding {

	namespace util {
		gsl_matrix_ptr generator_to_parity(const gsl_matrix &generator) {
		}

		gsl_matrix_ptr parity_to_generator(const gsl_matrix &parity) {
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

    void linear_code::sanity_check_properties() const {}

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
		: m_code{util::generator_to_standard_form(std::move(generator))}
		, m_basis_size{m_code->size1}
		, m_word_length{m_code->size2}
	{
		// NB: m_code is initialized first because it is declared first in the header file and not
		// 	   because it is put first in the list here.
	}

    linear_code::~linear_code() noexcept {
        gsl_matrix_free(m_code.get());
        if (m_decoding_table)
            gsl_matrix_free(m_decoding_table->get());
    }

    /*
     * Operations
     */

    gsl_vector_ptr linear_code::encode(const gsl_vector &word) {
        return with_redundancy(word);
    }

	gsl_vector_ptr linear_code::decode(const gsl_vector &) {}
}
