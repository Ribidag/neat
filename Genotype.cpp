#include "Genotype.hpp"

#include <iostream>

double Genotype::getRandomWeight() {
	return NumberGenerator::getSym(2.0);
}
double Genotype::getRandomBias() {
	return NumberGenerator::getSym(2.0);
}

Genotype::Genotype()
	: DirectedAcyclicGraph()
	, m_synapseGenes()
	, m_latestInnovation(0)
	, m_innovations()
	, m_neuronGenes() {}

Genotype::Genotype(const uint16_t numInputs, const uint16_t numOutputs)
	: DirectedAcyclicGraph()
	, m_synapseGenes()
	, m_latestInnovation(numInputs * numOutputs)
	, m_innovations()
	, m_neuronGenes() {

	m_innovations.reserve(numInputs * numOutputs);

	innov_t innovationNumber = 1;
	for (id_t inputId(-1); inputId >= -numInputs; --inputId) {
		m_neuronGenes.insert(std::make_pair(inputId, 0.0));

		for (id_t outputId(0); outputId < numOutputs; ++outputId) {
			addConnection(inputId, outputId);

			m_synapseGenes[innovationNumber] = SynapseGene(getRandomWeight(), inputId, outputId);
			m_innovations.push_back(innovationNumber);
			++innovationNumber;

			m_neuronGenes.insert(std::make_pair(outputId, getRandomBias()));
		}
	}

	orderNodes();
}

void Genotype::inheritSynapseGene(const innov_t innovationNumber, SynapseGene synapseGene) {
	const double enableRoll = NumberGenerator::getASym(1.0);
	if (enableRoll <= 0.25) {
		synapseGene.enable();
	}

	m_synapseGenes[innovationNumber] = synapseGene;
	m_innovations.push_back(innovationNumber);
	m_latestInnovation = std::max(m_latestInnovation, innovationNumber);

	const auto& inputOutputIds = synapseGene.getInputOutputIds();
	const id_t inputId = inputOutputIds.first;
	const id_t outputId = inputOutputIds.second;

	addConnection(inputId, outputId);
}

void Genotype::inheritNeuronGene(const id_t id, const double bias) {
	m_neuronGenes.insert(std::make_pair(id, bias));
}

Genotype::Genotype(const Genotype &fitterGenotype, const Genotype &weakerGenotype)
	: DirectedAcyclicGraph()
	, m_latestInnovation(0)
	, m_innovations() {

	const auto &fitterGenes = fitterGenotype.m_synapseGenes;
	const auto &weakerGenes = weakerGenotype.m_synapseGenes;
	const innov_t latestInnovationFitter = fitterGenotype.getLatestInnovation();
	const innov_t latestInnovationWeaker = weakerGenotype.getLatestInnovation();

	/* Derived and basal
	 *
	 * The greater the innovation number of a synapse gene, the more
	 * derived it is said to be. The genotype with the most derived
	 * synapse gene is the derived genotype. The other is the basal
	 * genotype. Note that the derived genotype is not necessarily
	 * the fitter genotype.
	 *
	 */
	innov_t maxInnovation, edgeInnovation;
	std::unordered_map<innov_t, SynapseGene> derivedGenes, basalGenes;
	bool derivedIsFitter;
	if (latestInnovationFitter > latestInnovationWeaker) {
		maxInnovation = latestInnovationFitter;
		edgeInnovation = latestInnovationWeaker;

		derivedGenes = fitterGenes;
		basalGenes = weakerGenes;

		derivedIsFitter = true;
	} else {
		maxInnovation = latestInnovationWeaker;
		edgeInnovation = latestInnovationFitter;

		derivedGenes = weakerGenes;
		basalGenes = fitterGenes;

		derivedIsFitter = false;
	}

	/* Inherit topology (synapse genes)
	 *
	 *
	 *
	 */
	for (innov_t i(0); i <= maxInnovation; ++i) {
		const auto& derivedIt = derivedGenes.find(i);
		const auto& basalIt = basalGenes.find(i);

		/* Excess Genes
		 *
		 * are only inherited if they come from the fitter parent.
		 *
		 */
		if (i > edgeInnovation && derivedIt != derivedGenes.end() && derivedIsFitter) {
			inheritSynapseGene(derivedIt->first, derivedIt->second);
			continue;
		}

		/* Matching Genes
		 *
		 * are picked randomly 50/50 from either parent.
		 * Therefore, it doesn't matter which is the fitter.
		 *
		 */

		if (derivedIt != derivedGenes.end() && basalIt != basalGenes.end()) {
			const double result = NumberGenerator::getASym(1.0);
			if (result < 0.5) {
				inheritSynapseGene(derivedIt->first, derivedIt->second);
			} else {
				inheritSynapseGene(basalIt->first, basalIt->second);
			}
			continue;
		}

		/* Disjoint Genes
		 *
		 * are only inherited if they come from the fitter parent.
		 *
		 */
		if (derivedIt != derivedGenes.end() && derivedIsFitter) {
			inheritSynapseGene(derivedIt->first, derivedIt->second);
		} else if (basalIt != basalGenes.end() && !derivedIsFitter) {
			inheritSynapseGene(basalIt->first, basalIt->second);
		}
	}

	/* Node order
	 *
	 * Since the structure of the new network is the same as that of the fitter parent,
	 * its node order can be transferred unchanged.
	 *
	 */
	setNodeOrder(fitterGenotype.getNodeOrder());

	/* Inherit neuron genes
	 *
	 * The biases are just transferred from the fitter parent since the structure,
	 * and thus the nodes, are the same.
	 *
	 */
	for (const id_t nodeId : getNodeOrder()) {
		const double bias = fitterGenotype.getNeuronGenes().find(nodeId)->second;
		inheritNeuronGene(nodeId, bias);
	}
}

void Genotype::addNeuronGene(const id_t id, const double bias) {
	if (!m_neuronGenes.count(id)) {
		m_neuronGenes.insert(std::make_pair(id, bias));
	} else {
		std::cerr << "Neuron gene " << id << " already exists! Cannot add it!" << std::endl;
	}
}

void Genotype::addSynapseGene(const innov_t innovationNumber, const double weight, const id_t inputId, const id_t outputId) {
	m_synapseGenes[innovationNumber] = SynapseGene(weight, inputId, outputId);
	m_innovations.push_back(innovationNumber);
	m_latestInnovation = std::max(m_latestInnovation, innovationNumber);

	if (!m_neuronGenes.count(inputId) || !m_neuronGenes.count(outputId)) {
		std::cerr << "Cannot add synapse gene between non-existent neuron genes" << std::endl;
	}

	addConnection(inputId, outputId);
}

void Genotype::splitSynapse(const innov_t firstInnovationNumber, const innov_t secondInnovationNumber, const id_t createdNeuronId, const innov_t synapseGeneToSplitId) {
	auto& synapseGeneToSplit = m_synapseGenes.find(synapseGeneToSplitId)->second;
	const std::pair<id_t, id_t> inputOutputIds = synapseGeneToSplit.getInputOutputIds();
	const id_t inputId = inputOutputIds.first;
	const id_t outputId = inputOutputIds.second;

	addNeuronGene(createdNeuronId, getRandomBias());

	addSynapseGene(firstInnovationNumber, synapseGeneToSplit.getWeight(), inputId, createdNeuronId);
	addSynapseGene(secondInnovationNumber, getRandomWeight(), createdNeuronId, outputId);

	synapseGeneToSplit.disable();
}

void Genotype::setSynapseWeight(const innov_t innovationNumber, const double weight) {
	m_synapseGenes[innovationNumber].setWeight(weight);
}

void Genotype::setNeuronBias(const id_t neuronGeneId, const double bias) {
	m_neuronGenes[neuronGeneId] = bias;
}

innov_t Genotype::findSplittableSynapse() const {
	const innov_t randomInnovation = m_innovations[std::floor(NumberGenerator::getASym(m_innovations.size()))];
	return randomInnovation;
}

const std::pair<id_t, id_t> Genotype::findConnectableNeurons(const uint64_t numInputs) const {
	uint16_t attempts = 0;

	while (attempts < m_neuronGenes.size() * 2) {
		++attempts;
		/*
		 * The end node must not be an input. Since inputs are always at depth 0,
		 * they go at the beginning of the node order. Outputs however, can be at arbitrary
		 * depth and cannot be sidestepped so easily. Thus, the start node must be verified
		 * not to be an output.
		 */
		const uint64_t endNodeIndex = std::floor(NumberGenerator::getASym(getNodeOrder().size() - numInputs)) + numInputs;
		const uint64_t startNodeIndex = std::floor(NumberGenerator::getASym(getNodeOrder().size()));

		const id_t startNodeId = getNodeOrder()[startNodeIndex];
		const id_t endNodeId = getNodeOrder()[endNodeIndex];

	//	std::cout << "Attempting to grow synapse " << startNodeId << " -> " << endNodeId;

		if (getNodes().find(startNodeId)->second.outputNodeIds.size() == 0) {
	//		std::cout << ". " << startNodeId << " was an output" << std::endl;
			continue;
		}

		if (connectable(startNodeId, endNodeId)) {
	//		std::cout << ". Connectable!" << std::endl;
			return std::make_pair(startNodeId, endNodeId);
		} else {
	//		std::cout << ". Not connectable!" << std::endl;
		}
	}

	return std::make_pair(-numInputs-1, -numInputs-1);
}
