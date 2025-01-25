#ifndef INCLUDE_NUMBERGENERATOR_HPP_
#define INCLUDE_NUMBERGENERATOR_HPP_

#include <random>

class NumberGenerator {
private:
    static std::random_device rd;
    static std::mt19937 gen;
    static std::uniform_real_distribution<float> sym_distribution;
    static std::uniform_real_distribution<float> asym_distribution;

public:
    NumberGenerator() = delete; // Prevent instantiation

    static void initialize() {
        gen = std::mt19937(rd());
        sym_distribution = std::uniform_real_distribution<float>(-1.0, 1.0);
        asym_distribution = std::uniform_real_distribution<float>(0.0, 1.0);
    }

    static double getSym(double range) {
        return range * sym_distribution(gen);
    }

    static double getASym(double range) {
        return range * asym_distribution(gen);
    }
};

#endif /* INCLUDE_NUMBERGENERATOR_HPP_ */
