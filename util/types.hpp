#ifndef NEAT_UTIL_TYPES_HPP_
#define NEAT_UTIL_TYPES_HPP_

#include <cstdint>
#include <limits>
#include <utility>

using innov_t = uint32_t;
using id_t = int;

struct SzudzikHash {
    std::size_t operator()(const std::pair<uint32_t, uint32_t>& p) const {
        constexpr uint32_t N = std::numeric_limits<uint32_t>::max();
        return p.first + p.second * N;
    }
};


#endif /* NEAT_UTIL_TYPES_HPP_ */
