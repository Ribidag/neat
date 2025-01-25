#ifndef NEAT_PHENOTYPE_HPP_
#define NEAT_PHENOTYPE_HPP_

#include <memory>

#include "DirectedAcyclicGraph.hpp"

class Genotype;

#include <iostream>

class Phenotype {
private:
	struct ExecutableNode {
		const std::unordered_set<id_t>& outputNodeIds;
		double bias;
		std::unordered_map<id_t, double> weights;
		std::unordered_set<id_t> enabledConnections;

		double input;

		ExecutableNode(const DirectedAcyclicGraph::Node& node) : outputNodeIds(node.outputNodeIds), bias(0), weights(), enabledConnections(), input(0) {}
	};

	std::unordered_map<id_t, std::unique_ptr<ExecutableNode>> m_executableNodes;
	std::vector<id_t>& m_nodeOrder;

public:
	Phenotype(Genotype& genotype);
	Phenotype(Phenotype&& other) noexcept;
	Phenotype& operator=(Phenotype&& other) noexcept;

	const std::unordered_map<id_t, double> execute(const std::unordered_map<id_t, double> &inputs);

	static double activateSigmoid(const double input);
};

#endif /* NEAT_PHENOTYPE_HPP_ */
