#include <memory>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <sstream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#define FMT_HEADER_ONLY
#include <fmt/core.h>

#ifndef LCODES_GSL_WRAPPER_H
#define LCODES_GSL_WRAPPER_H

namespace gsl_wrapper {
    // Get the niceties of modern smart pointers in this C API.
    using gsl_matrix_ptr = std::unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)>;
    using gsl_vector_ptr = std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)>;

    template<const int N>
    gsl_vector_ptr make_word(std::bitset<N> bits) {
        gsl_vector_ptr vec{gsl_vector_alloc(bits.size()), &gsl_vector_free};
        for (int i = 0; i < bits.size(); ++i)
            gsl_vector_set(vec.get(), i, bits.test(i));
        return vec;
    }

    void make_word(std::uint64_t bits, gsl_vector &vec);

    gsl_matrix_ptr gsl_matrix_transpose_non_square(const gsl_matrix &mx);

    gsl_vector_ptr gsl_matrix_mul_vector(const gsl_vector &vec, const gsl_matrix &matrix);

    template<typename... Args>
    void gsl_matrix_set_all(gsl_matrix &mx, Args... args) {
        const double arr[] = {static_cast<double>(args)...};
        size_t idx = 0;
        for (size_t i = 0; i < mx.size1; i++) {
            for (size_t j = 0; j < mx.size2; j++) {
                gsl_matrix_set(&mx, i, j, arr[idx++]);
            }
        }
    }

}

template<>
struct fmt::formatter<gsl_matrix>: formatter<string_view> {
    template<typename FormatContext>
    auto format(const gsl_matrix& m, FormatContext& ctx) const {
        std::stringstream sstr;
        for (size_t i = 0; i < m.size1; i++) {
            if (i > 0)
                sstr << '\n';
            for (size_t j = 0; j < m.size2; j++) {
                if (j > 0)
                    sstr << ' ';
                sstr << gsl_matrix_get(&m, i, j);
            }
        }
        return formatter<string_view>::format(sstr.str(), ctx);
    }
};

template<>
struct fmt::formatter<gsl_vector>: formatter<string_view> {
    template<typename FormatContext>
    auto format(const gsl_vector& m, FormatContext& ctx) const {
        std::stringstream sstr;
        sstr << '(';
        for (size_t i = 0; i < m.size - 1; i++)
            sstr << gsl_vector_get(&m, i) << ", ";
        sstr << gsl_vector_get(&m, m.size - 1);
        return formatter<string_view>::format(sstr.str(), ctx);
    }
};

#endif // LCODES_GSL_WRAPPER_H
