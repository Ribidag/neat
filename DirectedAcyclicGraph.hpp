#ifndef NEAT_DIRECTEDACYCLICGRAPH_HPP_
#define NEAT_DIRECTEDACYCLICGRAPH_HPP_

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "util/types.hpp"

class DirectedAcyclicGraph {
public:
	struct Node {
		uint32_t numInputs;
		std::unordered_set<id_t> outputNodeIds;
		uint32_t depth;

		Node() : numInputs(0), outputNodeIds(), depth(0) {}
	};

private:
	std::unordered_map<id_t, Node> m_nodes;
	std::vector<id_t> m_nodeOrder;

private:
	bool isFirstParentOfSecond(id_t parentId, id_t childId) const;
	bool isFirstAncestorOfSecond(id_t ancestorId, id_t descendantId) const;
	void assignNodeDepths();

public:
	bool connectable(id_t startId, id_t endId) const;
	void addConnection(id_t startId, id_t endId);

	void removeConnection(id_t startId, id_t endId); // unnecessary
	void splitConnection(const id_t startId, const id_t endId, const id_t id); // unnecessary

	void orderNodes();

	const std::vector<id_t>& getNodeOrder() const {
		return m_nodeOrder;
	}

	std::vector<id_t>& getNodeOrder() {
		return m_nodeOrder;
	}

	void setNodeOrder(const std::vector<id_t> nodeOrder) {
		m_nodeOrder = nodeOrder;
	}

	void addNode(const id_t id);

	const std::unordered_map<id_t, Node>& getNodes() const {
		return m_nodes;
	}

	const id_t getRandomNodeId() const;
};

#endif /* NEAT_DIRECTEDACYCLICGRAPH_HPP_ */
