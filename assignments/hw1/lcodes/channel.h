#ifndef CODINGTHEORY_HW1_CHANNEL_H
#define CODINGTHEORY_HW1_CHANNEL_H

#include <cstdint>

#include <gsl/gsl_vector.h>
#include <bitset>

#include "lcodes.h"

namespace coding {

    class channel {

    public: // Helpers
        template<const int N>
        static gsl_vector_ptr make_word(std::bitset<N> bits) {
            gsl_vector_ptr vec{gsl_vector_alloc(bits.size()), &gsl_vector_free};
            for (int i = 0; i < bits.size(); ++i)
                gsl_vector_set(vec.get(), i, bits.test(i));
            return vec;
        }

        static gsl_vector &with_noise(gsl_vector &word);

    public: // Instance control
        explicit channel(linear_code &&t_ecc)
                : m_ecc{std::move(t_ecc)} {}

    public: // Operations
        gsl_vector_ptr transfer(gsl_vector_ptr &&word);

    private:
        linear_code m_ecc;
    };

}

#endif // CODINGTHEORY_HW1_CHANNEL_H
