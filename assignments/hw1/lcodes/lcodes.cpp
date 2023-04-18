#include <string_view>
#include <algorithm>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>

#include "lcodes.h"

namespace coding {

    /*
     * Utilities
     */

    static void gsl_print_matrix(const gsl_matrix &matrix) {
        for (size_t i = 0; i < matrix.size1; i++) {
            for (size_t j = 0; j < matrix.size2; j++) {
                printf("%g", gsl_matrix_get(&matrix, i, j));
                if (j < matrix.size2 - 1)
                    printf(", ");
            }
            printf("\n");
        }
    }

    static gsl_matrix *gsl_copy_matrix(const gsl_matrix *proto) {
        auto *pmatrix = gsl_matrix_alloc(proto->size1, proto->size2);
        if (!pmatrix)
            throw linear_code_exception{"linear_code: gsl_matrix_alloc() returned nullptr"};
        if (gsl_matrix_memcpy(pmatrix, proto)) {
            gsl_matrix_free(pmatrix);
            throw linear_code_exception{"linear_code: gsl_matrix_memcpy() failed"};
        }
        return pmatrix;
    }

    struct hamming_metric {
    };
    using hamming = hamming_metric;

    template<typename T>
    std::size_t dist(bool, bool);

    template<typename T>
    std::size_t dist(const std::vector<bool> &, const std::vector<bool> &);

    template<>
    std::size_t dist<hamming>(bool a0, bool b0) { return a0 ^ b0; }

    template<>
    std::size_t dist<hamming>(const std::vector<bool> &a, const std::vector<bool> &b) {
        if (a.size() != b.size())
            throw linear_code_exception{"linear_code: Hamming distance called on words with different lengths"};

        std::size_t accumulated = 0;
        for (auto a_it = std::cbegin(a), b_it = std::cbegin(b); a_it != std::cend(a); ++a_it, ++b_it)
            accumulated += dist<hamming>(*a_it, *b_it);
        return accumulated;
    }

    static std::size_t weight(const std::vector<bool> &a) {
        return std::count_if(a.begin(), a.end(), [](bool b) { return b; });
    }

    static std::size_t weight(const gsl_vector_const_view &vec_cview) {
        std::size_t count_of_positives = 0;
        const auto &vec = vec_cview.vector;
        for (size_t i = 0; i < vec.size; ++i) {
            if (gsl_vector_get(&vec, i) >= std::numeric_limits<double>::epsilon())
                ++count_of_positives;
        }

        return count_of_positives;
    }

    /*
     * Helpers
     */

    void linear_code::evaluate_properties() {
        evaluate_min_distance();
        sanity_check_properties(); // Throws on error!

        // The capabilities property depends on the "core" ones. Therefore, we could postpone its
        // evaluation until after we have confirmed the validity of these other "core" properties.
        m_capabilities = (m_min_distance - 1) / 2;
    }

    void linear_code::evaluate_min_distance() {
        auto *reduced = gsl_copy_matrix(m_code);
        auto *permutation = gsl_permutation_alloc(m_basis_size);
        if (!permutation) {
            gsl_matrix_free(reduced);
            throw linear_code_exception{"linear_code: gsl_permutation_alloc() failed"};
        }
        int permutation_sign;

        if (gsl_linalg_LU_decomp(reduced, permutation, &permutation_sign) > 0) {
            gsl_matrix_free(reduced);
            throw linear_code_exception{"linear_code: gsl_linalg_LU_decomp() failed"};
        }

        auto min_distance_candidate = std::numeric_limits<std::size_t>::max();
        for (size_t i = 0; i < reduced->size1; i++) {
            gsl_vector_const_view row = gsl_matrix_const_row(reduced, i);
            if (const auto wt = weight(row); wt <= min_distance_candidate)
                min_distance_candidate = wt;
        }

        gsl_matrix_free(reduced);

        // Ta-daah!
        m_min_distance = min_distance_candidate;
    }

    void linear_code::evaluate_decoding_table() {}

    void linear_code::sanity_check_properties() const {
        static constexpr std::size_t MAX_WORD_LENGTH = 1 << 8;
        static constexpr std::size_t MAX_BASIS_SIZE = 1 << 8;
        const bool flag = m_basis_size <= MAX_BASIS_SIZE
                          && m_word_length <= MAX_WORD_LENGTH
                          && m_min_distance >= 1;

        if (!flag)
            throw linear_code_exception{"linear_code: sanity_check_properties() failed."};
    }

    /*
     * Instance control
     */

    linear_code::linear_code(gsl_matrix *t_code, size_t t_word_length, size_t t_basis_size)
            : m_code{t_code}, m_basis_size{t_basis_size}, m_word_length{t_word_length} {
        evaluate_properties();
    }

    linear_code::linear_code(const linear_code &rhs)
            : m_capabilities{rhs.m_capabilities}, m_min_distance{rhs.m_min_distance}, m_basis_size{rhs.m_basis_size},
              m_word_length{rhs.m_word_length} {
        m_code = gsl_copy_matrix(rhs.m_code);
        if (rhs.m_decoding_table)
            m_decoding_table = gsl_copy_matrix(rhs.m_decoding_table);
    }

    linear_code &linear_code::operator=(const linear_code &rhs) {
        m_capabilities = rhs.m_capabilities;
        m_min_distance = rhs.m_min_distance;
        m_basis_size = rhs.m_basis_size;
        m_word_length = rhs.m_word_length;
        if (m_code) {
            gsl_matrix_free(m_code);
            m_code = gsl_copy_matrix(rhs.m_code);
        }
        if (m_decoding_table) {
            gsl_matrix_free(m_decoding_table);
            m_decoding_table = gsl_copy_matrix(rhs.m_decoding_table);
        }
        return *this;
    }

    class linear_code
    linear_code::from_generator_matrix_as_array(const std::vector<double> &t_matrix, std::size_t t_word_size,
                                                std::size_t t_rows) {
        const std::size_t t_cols = t_matrix.size() / t_rows;
        gsl_matrix_const_view generator_matrix_const_view = gsl_matrix_const_view_array(t_matrix.data(), t_rows,
                                                                                        t_cols);
        auto *pmatrix = gsl_copy_matrix(&generator_matrix_const_view.matrix);

        return linear_code{pmatrix, t_word_size, t_rows};
    }

    class linear_code linear_code::from_dual(const linear_code &) {
    }

    linear_code::~linear_code() noexcept {
        gsl_matrix_free(m_code);
        if (m_decoding_table)
            gsl_matrix_free(m_decoding_table);
    }

    /*
     * Operations
     */

    linear_code::code_vector linear_code::encode(const linear_code::data_vector &) const noexcept {}

    linear_code::data_vector linear_code::decode(const linear_code::code_vector &) const {}
}