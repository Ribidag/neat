#include "SynapseGene.hpp"

SynapseGene::SynapseGene() : m_weight(0.0), m_enabled(true), m_inputOutputIds(std::make_pair(0, 0)) {}

SynapseGene::SynapseGene(const double weight, const id_t inputId, const id_t outputId)
	: m_weight(weight)
	, m_enabled(true)
	, m_inputOutputIds(std::make_pair(inputId, outputId)) {}
