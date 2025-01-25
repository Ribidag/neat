#ifndef NEAT_SYNAPSEGENE_HPP_
#define NEAT_SYNAPSEGENE_HPP_

#include <utility>

#include "util/types.hpp"

class SynapseGene {
private:
	double m_weight;
	bool m_enabled;

	std::pair<id_t, id_t> m_inputOutputIds;

public:
	SynapseGene();
	SynapseGene(const double weight, const id_t inputId, const id_t outputId);

	double getWeight() const {
		return m_weight;
	}

	void setWeight(double weight) {
		m_weight = weight;
	}

	void enable() {
		m_enabled = true;
	}

	void disable() {
		m_enabled = false;
	}

	const bool isEnabled() const {
		return m_enabled;
	}

	const std::pair<id_t, id_t>& getInputOutputIds() const {
		return m_inputOutputIds;
	}
};

#endif /* NEAT_SYNAPSEGENE_HPP_ */
