#include <array>
#include <random>

#include <gsl/gsl_matrix.h>

#include "lcodes.h"
#include "channel.h"

namespace gsl_wrapper {
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

int main() {
    using namespace coding;

    auto hamming_generator = []() {
        coding::gsl_matrix_ptr generator{gsl_matrix_alloc(4, 7), &gsl_matrix_free};
        gsl_wrapper::gsl_matrix_set_all(*generator.get(),
                                        1, 0, 0, 0, 0, 1, 1,
                                        0, 1, 0, 0, 1, 0, 1,
                                        0, 0, 1, 0, 1, 1, 0,
                                        0, 0, 0, 1, 1, 1, 1);
        return generator;
    }();
    auto hamming_code = linear_code::from_generator(std::move(hamming_generator));
    channel channel{std::move(hamming_code)};

    channel.transfer(channel::make_word<4>(0b1000));

    return 0;
}
