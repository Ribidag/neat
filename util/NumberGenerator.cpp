#include "NumberGenerator.hpp"

std::random_device NumberGenerator::rd;
std::mt19937 NumberGenerator::gen;
std::uniform_real_distribution<float> NumberGenerator::sym_distribution;
std::uniform_real_distribution<float> NumberGenerator::asym_distribution;
