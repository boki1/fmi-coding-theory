#include <iostream>
#include <array>
#include <cstddef>
#include <vector>
#include <random>

#include "lcodes.h"

using namespace coding;

namespace demo {

    template<const int N>
    std::vector<bool> make_word(uint64_t big_val) {
        std::vector<bool> word;
        word.reserve(N);
        for (int i = 0; i < N; ++i)
            word.push_back(big_val & (1 << i));
        return word;
    }

    template<typename T>
    // requires coding::is_correcting_code<T>
    class noisy_channel {

    public: // Instance control
        explicit noisy_channel(T *t_correcting_code)
                : m_correcting_code{t_correcting_code} {}

    public: // Operations
        void damage(std::vector<bool> &word) {
            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<> dist_num_errors(0, m_correcting_code->capabilities() - 1);
            std::vector<bool> flipped(word.size(), false);

            const auto errors = dist_num_errors(mt);
            std::uniform_int_distribution<int> dist_indices(0, word.size());

            for (int i = 0; i < errors; ++i) {
                int index = dist_indices(mt);
                while (flipped[index])
                    index = dist_indices(mt);
                word[index] = !word[index];
                flipped[index] = true;
            }
        }

        std::vector<bool> transfer(const std::vector<bool> &word) {
            auto tx_word = m_correcting_code->encode(word);
            damage(tx_word);
            const auto rx_word = m_correcting_code->decode(tx_word);
            return rx_word;
        }

    private:
        T *m_correcting_code;
    };

    void test_linear_5_2() {
        std::vector<double> generator_raw{
#include "generators/generator_matrix_for_5-2_lcode.inc"
        };

        try {
            auto lcode = linear_code::from_generator_matrix_as_array(generator_raw, 5, 2);

            noisy_channel<linear_code> channel{&lcode};
            auto word1 = make_word<5>(0b11001);
            channel.damage(word1);
        } catch (coding::linear_code_exception &e) {
            std::cout << e.what() << '\n';
        }
    }

    // The [8,4,3] Hamming code, which encodes 4 information bits into 8 bits by adding 4 parity bits.
    void test_linear_8_4() {
        std::vector<double> generator{
#include "generators/generator_matrix_for_8-4_lcode.inc"
        };

        try {
            auto lcode = linear_code::from_generator_matrix_as_array(generator, 8, 4);
        } catch (coding::linear_code_exception &e) {
            std::cout << e.what() << '\n';
        }
    }
}

int main() {
    using namespace demo;

    test_linear_5_2();
    test_linear_8_4();

    return 0;
}
