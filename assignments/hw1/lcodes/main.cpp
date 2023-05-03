#include <array>

#include <gsl/gsl_matrix.h>

#include "lcodes.h"
#include "channel.h"

static void separator(std::string_view prefix) {
    static constexpr std::string_view sep{"------------"};
    fmt::print("\n{}+ {} +{}\n", sep, prefix, sep);
}

namespace coding {
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
}

namespace coding {
    namespace test {
        void hamming() {
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
            hamming_code.set_name("Hamming [7, 4]");

            fmt::print("Generator for Hamming [7, 4]:\n{}\n", hamming_code.code());

            auto hamming_parity = util::generator_to_parity(hamming_code.code());
            fmt::print("Parity-check matrix:\n{}\n", *hamming_parity);

            auto hamming_generator_ = util::parity_to_generator(*hamming_parity);
            fmt::print("Generator from parity matrix:\n{}\n", *hamming_generator_);

            fmt::print("{}", hamming_code);

//    channel channel{std::move(hamming_code)};
//
//    channel.transfer(channel::make_word<4>(0b1000));
        }

        void lin_4_2() {
            auto lin_4_2_generator = []() {
                coding::gsl_matrix_ptr generator{gsl_matrix_alloc(2, 4), &gsl_matrix_free};
                gsl_wrapper::gsl_matrix_set_all(*generator.get(),
                                                1, 0, 1, 1,
                                                0, 1, 0, 1);
                return generator;
            }();

            auto lin_4_2_code = linear_code::from_generator(std::move(lin_4_2_generator));
            fmt::print("Generator for C[4, 2]:\n{}\n", lin_4_2_code.code());

            auto lin_4_2_parity = util::generator_to_parity(lin_4_2_code.code());
            fmt::print("Parity-check matrix:\n{}\n", *lin_4_2_parity);

            auto lin_4_2_generator_ = util::parity_to_generator(*lin_4_2_parity);
            fmt::print("Generator from parity matrix:\n{}\n", *lin_4_2_generator_);
        }

        void lin_6_3() {
            auto lin_6_3_generator = []() {
                coding::gsl_matrix_ptr generator{gsl_matrix_alloc(3, 6), &gsl_matrix_free};
                gsl_wrapper::gsl_matrix_set_all(*generator.get(),
                                                1, 0, 0, 0, 1, 1,
                                                0, 1, 0, 1, 0, 1,
                                                0, 0, 1, 1, 1, 0);
                return generator;
            }();

            auto lin_6_3_code = linear_code::from_generator(std::move(lin_6_3_generator));
            fmt::print("Generator for C[6, 3]:\n {}\n", lin_6_3_code.code());

            auto lin_6_3_parity = util::generator_to_parity(lin_6_3_code.code());
            fmt::print("Parity-check matrix:\n {}\n", *lin_6_3_parity);

            auto lin_6_3_generator_ = util::parity_to_generator(*lin_6_3_parity);
            fmt::print("Generator from parity matrix:\n{}\n", *lin_6_3_generator_);
        }

        void lin_5_3() {
            auto lin_5_3_generator = []() {
                coding::gsl_matrix_ptr generator{gsl_matrix_alloc(3, 5), &gsl_matrix_free};
                gsl_wrapper::gsl_matrix_set_all(*generator.get(),
                                                1, 0, 0, 1, 1,
                                                0, 1, 0, 0, 1,
                                                0, 0, 1, 1, 1);
                return generator;
            }();

            auto lin_5_3_code = linear_code::from_generator(std::move(lin_5_3_generator));
            fmt::print("Generator for C[5, 3]:\n{}\n", lin_5_3_code.code());

            auto lin_5_3_parity = util::generator_to_parity(lin_5_3_code.code());
            fmt::print("Parity-check matrix:\n{}\n", *lin_5_3_parity);

            auto lin_5_3_generator_ = util::parity_to_generator(*lin_5_3_parity);
            fmt::print("Generator from parity matrix:\n{}\n", *lin_5_3_generator_);
        }
    }
}

int main() {
    using namespace coding::test;

    try {

    separator("Hamming code [7, 4]");
    hamming();
//
//    separator("Linear code [4, 2]");
//    lin_4_2();

//    separator("Linear code [6, 3]");
//    lin_6_3();

//    separator("Linear code [5, 3]");
//    lin_5_3();

    } catch (coding::linear_code_exception linear_code_exception) {
        fmt::print("{}", linear_code_exception.what());
        return 1;
    }

    return 0;
}
