#include <cstdint>
#include <optional>

#ifndef CODINGTHEORY_HW1_LCODES_H
#define CODINGTHEORY_HW1_LCODES_H

namespace coding {

/// Implementation of binary linear code.
/// TODO: Describe better :).
template <uint16_t N, uint16_t K>
class linear_code {
public:
    using code_vector = uint32_t;
    using data_vector = uint32_t;

    [[nodiscard]]
    code_vector encode(const data_vector&) const noexcept { /* FIXME: Implement */ }

    [[nodiscard]]
    std::optional<data_vector> decode(const code_vector &) const noexcept { /* FIXME: Implement */ }
public:
    // FIXME:
    using generator_matrix = Eigen::MatrixXd;
    using check_matrix = Eigen::MatrixXd;

    static linear_code<N, K> from_generator_matrix(const generator_matrix &) { /* FIXME: Implement */ }
    static linear_code<N, K> from_check_matrix(const check_matrix &) { /* FIXME: Implement */ }
    static linear_code<N, K> from_dual(const linear_code<N, K> &) { /* FIXME: Implement */ }

public: // Basic properties of the code.
    [[nodiscard]]
    std::uint32_t word_length() const noexcept { return N; }

    [[nodiscard]]
    std::uint32_t basis_size() const noexcept { return K; }

    [[nodiscard]]
    std::uint32_t num_codes() const noexcept { return 1 << K; }

    [[nodiscard]]
    std::uint32_t min_distance() const noexcept { /* FIXME: Implement */ }

    [[nodiscard]]
    std::uint32_t capabilities() const noexcept { /* FIXME: Implement */ }

    [[nodiscard]]
    std::uint32_t radius() const noexcept { /* FIXME: Implement */ }
};

/// Emulates a channel with noise.
class channel {

};

}

#endif // CODINGTHEORY_HW1_LCODES_H
