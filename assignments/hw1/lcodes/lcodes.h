#include <cstdint>
#include <optional>
#include <utility>
#include <exception>
#include <string>
#include <vector>

#include <gsl/gsl_matrix.h>

#ifndef CODINGTHEORY_HW1_LCODES_H
#define CODINGTHEORY_HW1_LCODES_H

namespace coding {

template<typename T>
concept is_correcting_code = requires(const T& code, const std::vector<bool>& message) {
    { code.encode(message) } -> std::same_as<const std::vector<bool>>;
    { code.decode(message) } -> std::same_as<const std::vector<bool>>;
};


class linear_code_exception : std::exception {
    std::string m_msg;
public:
    explicit linear_code_exception(std::string msg) : m_msg{std::move(msg)} {}

    [[nodiscard]]
    const char *what() const noexcept override { return m_msg.c_str(); }
};

/// Implementation of binary linear code.
/// TODO: Describe better :).
class linear_code {
    // Caches the properties of the code in the corresponding member variables.
    // This is called once right after construction.
    void evaluate_properties();

    void evaluate_min_distance();

    void evaluate_decoding_table();

    // Throws on invalid data.
    void sanity_check_properties() const;

public: // Operations:
    using code_vector = std::vector<bool>;
    using data_vector = std::vector<bool>;

    // Encodes a single data vector into a code vector.
    [[nodiscard]] linear_code::code_vector encode(const data_vector&) const noexcept;

    // Tries to decode a single code vector into a data vector.
    // On failure, a linear_code_exception is thrown.
    [[nodiscard]] linear_code::data_vector decode(const code_vector &) const;

public: // Instance control:
    // TODO: Make movable.

    explicit linear_code(gsl_matrix *t_code, size_t t_word_length, size_t t_basis_size);
    linear_code(const linear_code &);
    linear_code& operator=(const linear_code&);

    ~linear_code() noexcept;

    // Rationale for this ugly parameter list:
    // Although `lcodes` depends on GSL, I didn't want any client code to have to construct any GSL objects in order
    // to be able to use it. The goal here would be to not force `lcodes-demo` to link explicitly with GSl. Also,
    // the simplistic approach isn't thaaaat terrible.

    // FIXME: We do not use double. Improve upon the memory layout by using `gsl_matrix_uint`.

    [[nodiscard]]
    static linear_code from_generator_matrix_as_array(const std::vector<double> &, std::size_t, std::size_t);

    [[nodiscard]]
    static linear_code from_dual(const linear_code &);

public: // Basic properties of the code:
    [[nodiscard]]
    std::size_t word_length() const noexcept { return m_word_length; }

    [[nodiscard]]
    std::size_t basis_size() const noexcept { return m_basis_size; }

    [[nodiscard]]
    std::size_t num_codes() const noexcept { return 1 << m_basis_size; }

    [[nodiscard]]
    std::size_t min_distance() const noexcept { return m_min_distance; }

    [[nodiscard]]
    std::size_t capabilities() const noexcept { return m_capabilities; }

private:
    // m_word_length is the size of each code word.
    std::size_t m_word_length;

    // m_basis_size is the number of elements in the basis of the vector subspace that is the code.
    std::size_t m_basis_size;

    // m_min_distance in the minimal distance between two code words (using the Hamming metric).
    // It is defined as d(C) = min { d(a, b) | a, b \in C }. 0 is obviously an invalid value.
    std::size_t m_min_distance{0};

    // m_capabilities stores the unambiguous decoding capabilities of the code, i.e - the maximum
    // amount of errors in a single word that code may correct properly.
    std::size_t m_capabilities{0};

    // m_code is a matrix with elements of GF(2) with dimensions m_basis_size x (m_basis_size + m_word_length).
    gsl_matrix *m_code{nullptr};

    // m_decoding_table is what is known as Slepian's table and is one of the methods for decoding linear codes.
    // The member gets filled once at the first time of word decoding and gets used without modifications after that.
    gsl_matrix *m_decoding_table{nullptr};
};

}

#endif // CODINGTHEORY_HW1_LCODES_H
