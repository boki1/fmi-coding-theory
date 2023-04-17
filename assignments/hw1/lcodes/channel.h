#include <cstddef>

#ifndef CODINGTHEORY_HW1_CHANNEL_H
#define CODINGTHEORY_HW1_CHANNEL_H

namespace coding {

// Channel with noise.
class channel {
    std::byte transfer(const std::byte &);
};

}

#endif // CODINGTHEORY_HW1_CHANNEL_H
