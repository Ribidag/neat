#ifndef NEAT_UTIL_RANDOMUNORDEREDSET_HPP_
#define NEAT_UTIL_RANDOMUNORDEREDSET_HPP_

#include <vector>

#include "../../snake_game/include/NumberGenerator.hpp"

template<typename T>
class RandomUnorderedSet {
private:
	std::vector<T> container;

public:
	RandomUnorderedSet() : container() {}
	RandomUnorderedSet(const size_t initialSize) : container(initialSize) {}

	void insert(const T value) {
		container.push_back(value);
	}
	void erase(const T value);
	void clear();
	void size();
};

#endif /* NEAT_UTIL_RANDOMUNORDEREDSET_HPP_ */
