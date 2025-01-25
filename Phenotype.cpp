#include "Phenotype.hpp"

#include <cmath>
#include <vector>

#include "Genotype.hpp"

#include <iostream>

Phenotype::Phenotype(Genotype &genotype)
	: m_executableNodes()
	, m_nodeOrder(genotype.getNodeOrder()) {

	const auto& nodes = genotype.getNodes();

	// Create executable nodes
	for (const auto& nodeIt : nodes) {
		const id_t nodeId = nodeIt.first;
		const auto& node = nodeIt.second;

		m_executableNodes[nodeId] = std::make_unique<ExecutableNode>(node);
	}

	// store synapse data
	for (const auto& synapseGeneIt : genotype.getSynapseGenes()) {
		const auto& synapseGene = synapseGeneIt.second;
		const id_t startNodeId = synapseGene.getInputOutputIds().first;
		const id_t endNodeId = synapseGene.getInputOutputIds().second;
		const double weight = synapseGene.getWeight();
		const bool enabled = synapseGene.isEnabled();

		auto& startExecutableNode = m_executableNodes[startNodeId];
		startExecutableNode->weights[endNodeId] = weight;

		if (enabled) {
			startExecutableNode->enabledConnections.insert(endNodeId);
		}
	}

	// store neuron data
	for (const auto& neuronGeneIt : genotype.getNeuronGenes()) {
		const id_t neuronGeneId = neuronGeneIt.first;
		const double bias = neuronGeneIt.second;

		m_executableNodes[neuronGeneId]->bias = bias;
	}
}

Phenotype::Phenotype(Phenotype&& other) noexcept
    : m_executableNodes(std::move(other.m_executableNodes)),
      m_nodeOrder(other.m_nodeOrder) {
}

Phenotype& Phenotype::operator=(Phenotype&& other) noexcept {
    if (this != &other) {
        m_executableNodes = std::move(other.m_executableNodes);
        m_nodeOrder = other.m_nodeOrder;
    }
    return *this;
}

double Phenotype::activateSigmoid(const double input) {
	return 1.0 / (1 + std::exp(-4.9 * input));
}

const std::unordered_map<id_t, double> Phenotype::execute(const std::unordered_map<id_t, double> &inputs) {
	//reset
	for (const id_t nodeId : m_nodeOrder) {
		const auto& executableNode = m_executableNodes[nodeId];
		executableNode->input = 0;
	}

	// load inputs
	for (const auto& inputIt : inputs) {
		const id_t inputNodeId = inputIt.first;
		const double inputValue = inputIt.second;

		m_executableNodes[inputNodeId]->input = inputValue;
//		std::cout << "Loaded " << m_executableNodes[inputNodeId]->input << " into node " << inputNodeId << std::endl;
	}

	//execute
	std::unordered_map<id_t, double> outputs;
	for (const id_t nodeId : m_nodeOrder) {
		const auto& executableNode = m_executableNodes[nodeId];

		const auto& outputNodeIds = executableNode->outputNodeIds;
		const auto& weights = executableNode->weights;
		const double rawOutput = executableNode->input + executableNode->bias;
		const auto& enabledConnections = executableNode->enabledConnections;

		const double output = (inputs.count(nodeId)) ? rawOutput : activateSigmoid(rawOutput);

//		std::cout << "Activated output from node " << nodeId << " is " << output << " using bias " << executableNode->bias << std::endl;

		if (outputNodeIds.size() == 0) {
//			std::cout << "Node " << nodeId << " is outputting " << rawOutput << std::endl;
			outputs[nodeId] = rawOutput;
			continue;
		}

		for (const id_t outputNodeId : outputNodeIds) {
			const double weight = weights.find(outputNodeId)->second;
			const bool enabled = enabledConnections.count(outputNodeId);

			if (enabled) {
				auto& outputExecutableNode = m_executableNodes[outputNodeId];
				const double input = weight * output;
				outputExecutableNode->input += input;

//				std::cout << "	Output node " << outputNodeId << " now has a total input of " << outputExecutableNode->input << " after an addition of " << input << " through connection weight " << weight << std::endl;
			} else {
//				std::cout << "	Connection to node " << outputNodeId << " was disabled!" << std::endl;
			}
		}
	}
//	std::cout << std::endl;

	return outputs;
}
