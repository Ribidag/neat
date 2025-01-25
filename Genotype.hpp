#ifndef NEAT_GENOTYPE_HPP_
#define NEAT_GENOTYPE_HPP_

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "DirectedAcyclicGraph.hpp"
#include "SynapseGene.hpp"
#include "util/NumberGenerator.hpp"

class Genotype : public DirectedAcyclicGraph {
private:
	std::unordered_map<innov_t, SynapseGene> m_synapseGenes;
	innov_t m_latestInnovation;
	std::vector<innov_t> m_innovations;

	std::unordered_map<id_t, double> m_neuronGenes;

private:
	void inheritSynapseGene(const innov_t innovationNumber, SynapseGene synapseGene);
	void inheritNeuronGene(const id_t id, const double bias);

public:
	Genotype();
	Genotype(const uint16_t numInputs, const uint16_t numOutputs);
	Genotype(const Genotype& fitterGenotype, const Genotype& weakerGenotype);

	void addNeuronGene(const id_t id, const double bias);
	void addSynapseGene(const innov_t innovationNumber, const double weight, const id_t inputId, const id_t outputId);

	void splitSynapse(const innov_t firstInnovationNumber, const innov_t secondInnovationNumber, const id_t createdNeuronId, const innov_t synapseGeneToSplitId);

	innov_t findSplittableSynapse() const;
	const std::pair<id_t, id_t> findConnectableNeurons(const uint64_t numInputs) const;

public:
	innov_t getLatestInnovation() const {
		return m_latestInnovation;
	}

	const std::vector<innov_t>& getInnovations() const {
		return m_innovations;
	}

	const std::unordered_map<innov_t, SynapseGene>& getSynapseGenes() const {
		return m_synapseGenes;
	}

	const std::unordered_map<id_t, double>& getNeuronGenes() const {
		return m_neuronGenes;
	}

	void setSynapseWeight(const innov_t innovationNumber, const double weight);
	void setNeuronBias(const id_t neuronGeneId, const double bias);

public:
	static double getRandomWeight();
	static double getRandomBias();
};

#endif /* NEAT_GENOTYPE_HPP_ */
