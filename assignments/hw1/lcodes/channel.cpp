#include <cassert>
#include <random>

#include "channel.h"

namespace coding {

    gsl_vector &channel::with_noise(gsl_vector &word) {
        auto toggle_bit_at = [&word](int index) {
            bool current = gsl_vector_get(&word, index) == 0;
            gsl_vector_set(&word, index, !current);
        };

        static std::random_device rd;
        static std::mt19937 mt(rd());
        std::uniform_int_distribution<> dist_num_errors(0, 10);
        std::vector<bool> flipped(word.size, false);

        const auto errors = dist_num_errors(mt);
        std::uniform_int_distribution<int> dist_indices(0, word.size);

        for (int i = 0; i < errors; ++i) {
            int index = dist_indices(mt);
            while (flipped[index])
                index = dist_indices(mt);
            toggle_bit_at(index);
            flipped[index] = true;
        }

        return word;
    }

    gsl_vector_ptr channel::transfer(coding::gsl_vector_ptr &&word) {
        auto tx_word = m_ecc.encode(*word);
        (void) with_noise(*tx_word);
        auto rx_word = m_ecc.decode(*tx_word);
        return rx_word;
    }
    
}
