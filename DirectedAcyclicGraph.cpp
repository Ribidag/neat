#include "DirectedAcyclicGraph.hpp"

#include <algorithm>
#include <stdexcept>

#include <iostream>

bool DirectedAcyclicGraph::isFirstParentOfSecond(id_t parentId, id_t childId) const {

	const std::unordered_set<id_t>& childIds = m_nodes.find(parentId)->second.outputNodeIds;
	return childIds.count(childId);
}

bool DirectedAcyclicGraph::isFirstAncestorOfSecond(id_t ancestorId, id_t descendantId) const {

	// base case
	if (isFirstParentOfSecond(ancestorId, descendantId)) {
		return true;
	}

	// recursion
	const std::unordered_set<id_t>& childIds = m_nodes.find(ancestorId)->second.outputNodeIds;
	for (const id_t childId : childIds) {
		if (isFirstAncestorOfSecond(childId, descendantId)) {
			return true;
		}
	}

	return false;
}

bool DirectedAcyclicGraph::connectable(id_t startId, id_t endId) const {


	if (!m_nodes.count(startId)) {
		return false;
	}

	if (!m_nodes.count(endId)) {
		return false;
	}

	if (startId == endId) {
		return false;
	}

	if (isFirstAncestorOfSecond(endId, startId)) {
		return false;
	}

	if (isFirstParentOfSecond(startId, endId)) {
		return false;
	}

	return true;
}

void DirectedAcyclicGraph::addConnection(id_t startId, id_t endId) {

	//std::cout << "	start: " << startId << ", end: " << endId << std::endl;

	addNode(startId);
	addNode(endId);

	//std::cout << "	num nodes: " << m_nodes.size() << std::endl;

	m_nodes.find(startId)->second.outputNodeIds.insert(endId);
	++m_nodes.find(endId)->second.numInputs;
}

void DirectedAcyclicGraph::removeConnection(id_t startId, id_t endId)
{
	auto& startingNode = m_nodes.find(startId)->second;
	if (startingNode.outputNodeIds.count(endId)) {
		startingNode.outputNodeIds.erase(endId);

		auto& endingNode = m_nodes.find(endId)->second;
		--endingNode.numInputs;
	}
}


void DirectedAcyclicGraph::splitConnection(const id_t startId, const id_t endId, const id_t id) {
	addConnection(startId, id);
	addConnection(id, endId);
}


void DirectedAcyclicGraph::assignNodeDepths() {

	std::vector<id_t> startNodeIds;

	std::unordered_map<id_t, uint32_t> numInputsAtId;

	for (auto &it : m_nodes) {
		auto& id = it.first;
		auto& n = it.second;

		numInputsAtId[id] = n.numInputs;

		if (n.numInputs == 0) {
			n.depth = 0;
			startNodeIds.push_back(id);
		}
	}

//	for (const auto& n : m_nodes) {
//		std::cout << "Node: " << n.first << " has " << n.second.numInputs << " inputs" << std::endl;
//	}
//
//	std::cout << "StartNodes:";
//	for (const auto& i : startNodeIds) {
//		std::cout << " " << i;
//	}
//	std::cout << std::endl;

	while (!startNodeIds.empty()) {
		const id_t startingNodeId = startNodeIds.back();
		startNodeIds.pop_back();
		const Node& startingNode = m_nodes.find(startingNodeId)->second;
		//std::cout << "Starting at node: " << m_nodes.find(startingNodeId)->first << " with " << startingNode.outputNodeIds.size() << " outputs!" << std::endl;

		for (const auto& outputNodeId : startingNode.outputNodeIds) {
			//std::cout << "	Outputing to node: " << outputNodeId<< std::endl;
			auto& outputNode = m_nodes.find(outputNodeId)->second;
			outputNode.depth = std::max(outputNode.depth, startingNode.depth + 1);
			--numInputsAtId.find(outputNodeId)->second;

			if (numInputsAtId.find(outputNodeId)->second == 0) {
				startNodeIds.push_back(outputNodeId);
			}
		}
	}

	//std::cout << "Assigned depths!" << std::endl << std::endl;
}

void DirectedAcyclicGraph::orderNodes() {
	assignNodeDepths();

	m_nodeOrder.clear();
	m_nodeOrder.reserve(m_nodes.size());

	for (const auto& it : m_nodes) {
		const id_t id = it.first;
		m_nodeOrder.push_back(id);
	}

	std::sort(m_nodeOrder.begin(), m_nodeOrder.end(), [this](id_t i, id_t j) {
		return m_nodes.find(i)->second.depth < m_nodes.find(j)->second.depth;
	});
}

void DirectedAcyclicGraph::addNode(const id_t id) {
	const auto& it = m_nodes.find(id);
	if (it == m_nodes.end()) {
		//std::cout << "		Created node: " << id << std::endl;
		m_nodes[id] = Node();
	}
}
