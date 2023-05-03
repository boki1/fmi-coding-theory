#include <cassert>
#include <cstdint>
#include <optional>
#include <utility>
#include <exception>
#include <string>
#include <vector>
#include <memory>
#include <sstream>
#include <bitset>

//! Sadly, there is a naming collision in between 2 of the 3 libraries that I am using.
//! GSL means both "Guidelines Support Library" and "GNU Scientific Library". Thus,
//! there are some places in the implementation where things just get weird. Eg: gsl::owner<gsl_matrix *> a;

#define FMT_HEADER_ONLY
#include <fmt/core.h>

// GNU Scientific Library
#include <gsl/gsl_matrix.h>

#ifndef CODINGTHEORY_HW1_LCODES_H
#define CODINGTHEORY_HW1_LCODES_H

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

namespace coding {

// Get the niceties of modern smart pointers in this C API.
    using gsl_matrix_ptr = std::unique_ptr<gsl_matrix, decltype(&gsl_matrix_free)>;
    using gsl_vector_ptr = std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)>;

    class linear_code_exception : std::exception {
        std::string m_msg;
    public:
        explicit linear_code_exception(std::string msg) : m_msg{std::move(msg)} {}

        [[nodiscard]]
        const char *what() const noexcept override { return m_msg.c_str(); }
    };

    namespace util {
        // Converts between different linear code representations - generator matrix and parity check equations.
        // Both functions return the resulting *newly allocated* matrix with the expected properties.
        gsl_matrix_ptr generator_to_parity(const gsl_matrix &);

        gsl_matrix_ptr parity_to_generator(const gsl_matrix &);

        // Given a generator matrix, puts it in standard form: G = (E|A).
        // Returns the same reference that was provided as an argument. Does not allocated anew.
        gsl_matrix &generator_to_standard_form(gsl_matrix &);
    }

    namespace gsl_wrapper {
        gsl_matrix_ptr gsl_matrix_transpose_non_square(const gsl_matrix &mx);

        template<const int N>
        static gsl_vector_ptr make_word(std::bitset<N> bits) {
            gsl_vector_ptr vec{gsl_vector_alloc(bits.size()), &gsl_vector_free};
            for (int i = 0; i < bits.size(); ++i)
                gsl_vector_set(vec.get(), i, bits.test(i));
            return vec;
        }

        // In-place dynamic variant of the other overload.
        static void make_word(std::uint64_t bits, gsl_vector &vec) {
            assert(vec.size <= sizeof(bits) * CHAR_BIT);
            for (int i = 0; i < vec.size; ++i)
                gsl_vector_set(&vec, i, bits & (1 << i));
        }
    }

/// Implementation of binary linear code.
/// TODO: Describe better :).
    class linear_code {

    private: // Helpers:
        // The following helper functions, cache the properties of the
        // code in the corresponding member variables.
        void evaluate_min_distance();

        void evaluate_capabilities();

        void evaluate_decoding_table();

        void sanity_check_properties() const; // Throws on invalid data.

        gsl_vector_ptr with_redundancy(const gsl_vector &word) const; // Throws on invalid input.

    public: // Operations:
        [[nodiscard]] static std::size_t weight(const gsl_vector &);

        // Encodes a single data vector into a code vector.
        [[nodiscard]] gsl_vector_ptr encode(const gsl_vector &);

        // Tries to decode a single code vector into a data vector.
        // On failure, a linear_code_exception is thrown.
        [[nodiscard]] gsl_vector_ptr decode(const gsl_vector &);

    private:
        // Construct codes only using the "named ctors".
        explicit linear_code(gsl_matrix_ptr &&);

    public: // Instance control:

        // Use only moves on instances of this class.
        linear_code(const linear_code &) = delete;

        linear_code &operator=(const linear_code &) = delete;

        linear_code(linear_code &&) noexcept = default;

        linear_code &operator=(linear_code &&) noexcept = default;

        [[nodiscard]] static linear_code from_parity_equations(gsl_matrix_ptr &&);

        [[nodiscard]] static linear_code from_generator(gsl_matrix_ptr &&);

        [[nodiscard]] static linear_code from_dual(const linear_code &);

    public: // Basic properties of the code:
        [[nodiscard]] std::size_t word_length() const noexcept { return m_word_length; }

        [[nodiscard]] std::size_t basis_size() const noexcept { return m_basis_size; }

        [[nodiscard]] std::size_t num_codes() const noexcept { return 1 << m_basis_size; }

        [[nodiscard]] std::size_t min_distance() const noexcept { return m_min_distance; }

        [[nodiscard]] std::size_t max_errors_detect() const noexcept { return m_max_errors_detect; }

        [[nodiscard]] std::size_t max_errors_correct() const noexcept { return m_max_errors_correct; }

        [[nodiscard]] std::size_t radius() const noexcept { return m_radius; }

        [[nodiscard]] const gsl_matrix &code() const noexcept { return *m_code; }

        void set_name(const std::string &name) noexcept { m_name = name; }

        [[nodiscard]] const std::string &name() const noexcept { return m_name; }

    private:
        std::string m_name;

        // m_code is a matrix with elements of GF(2) with dimensions m_basis_size x (m_basis_size + m_word_length).
        // We store the code in standard format generator matrix - G = (E|A).
        gsl_matrix_ptr m_code;

        // m_word_length is the size of each code word.
        std::size_t m_word_length;

        // m_basis_size is the number of elements in the basis of the vector subspace that is the code.
        std::size_t m_basis_size;

        // m_min_distance in the minimal distance between two code words (using the Hamming metric).
        // It is defined as d(C) = min { d(a, b) | a, b \in C }. 0 is obviously an invalid value.
        std::size_t m_min_distance{0};

        // m_max_errors_detect stores the unambiguous decoding capabilities of the code, i.e - the maximum
        // amount of errors in a single word that code may correct properly.
        std::size_t m_max_errors_detect{0};
        std::size_t m_max_errors_correct{0};

        std::size_t m_radius{0};

        // m_decoding_table is what is known as Slepian's table and is one of the methods for decoding linear codes.
        // The member gets filled once at the first time of word decoding and gets used without modifications after that.
        std::optional<gsl_matrix_ptr> m_decoding_table;
    };

}

template<>
struct fmt::formatter<coding::linear_code>: formatter<string_view> {
    template<typename FormatContext>
    auto format(const coding::linear_code& lcode, FormatContext& ctx) const {
        std::stringstream sstr;
        sstr << "linear_code <" << lcode.name() << "> {\n";
        sstr << "  - n: " << lcode.word_length() << '\n';
        sstr << "  - k: " << lcode.basis_size() << '\n';
        sstr << "  - d: " << lcode.min_distance() << '\n';
        sstr << "  - R: " << lcode.radius() << '\n';
        sstr << "  - max errors:\n  - detect: " << lcode.max_errors_detect() << "\n    - correct: " << lcode.max_errors_correct() << '\n';
        sstr << "};";
        return formatter<string_view>::format(sstr.str(), ctx);
    }
};

#endif // CODINGTHEORY_HW1_LCODES_H
