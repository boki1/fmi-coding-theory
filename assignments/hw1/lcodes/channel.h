#ifndef CODINGTHEORY_HW1_CHANNEL_H
#define CODINGTHEORY_HW1_CHANNEL_H

#include <cstdint>

#include <gsl/gsl_vector.h>
#include <bitset>

#include "lcodes.h"
#include "gsl_wrapper.h"

namespace coding {

    class channel {

    public: // Instance control
        explicit channel(linear_code &&t_ecc)
                : m_ecc{std::move(t_ecc)} {}

    public: // Operations
        gsl_wrapper::gsl_vector_ptr transfer(const gsl_vector &word);

        static gsl_vector &with_noise(gsl_vector &word);

    private:
        linear_code m_ecc;
    };

}

#endif // CODINGTHEORY_HW1_CHANNEL_H
