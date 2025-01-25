#include "DAGRenderer.hpp"

#include <unordered_map>
#include <utility>

#include "../DirectedAcyclicGraph.hpp"
#include "../util/NumberGenerator.hpp"

#include <iostream>

DAGRenderer::DAGRenderer()
	: m_position({0, 0})
	, m_nodeSpacing(0)
	, m_layerSpacing(0)
	, m_nodeRadius(0) {}

void DAGRenderer::addDrawableNode(const uint32_t depth, const id_t id) {
	const std::vector<DrawableNode>& nodesInLayer = m_nodes[depth];
	const uint32_t numNodesInLayer = nodesInLayer.size();

	double y = 0.0;
	if (numNodesInLayer > 0) {
	    if (numNodesInLayer & 1) { // odd
	        y = -m_nodeSpacing * (numNodesInLayer + 1) / 2;
	    } else {
	        y = m_nodeSpacing * numNodesInLayer / 2;
	    }
	}

	double depthDistance = depth*m_layerSpacing;

	const float randomHeightOffset = NumberGenerator::getSym(m_nodeSpacing*0.05);
	y += randomHeightOffset;
	const float randomDepthOffset = NumberGenerator::getSym(m_layerSpacing*0.1);
	depthDistance += randomDepthOffset;

	const sf::Vector2f nodePosition = m_position + sf::Vector2f(depthDistance, y);
	DrawableNode node(nodePosition, m_nodeRadius, id);

	m_nodes[depth].push_back(node);
}

void DAGRenderer::loadDAG(const DirectedAcyclicGraph &dag, const float width, const float height) {
	m_nodes.clear();
	m_connections.clear();

	m_position = sf::Vector2f(width * 0.15, height / 2.0);

	const auto &nodes = dag.getNodes();
	const auto &nodeOrder = dag.getNodeOrder();

	// Create layers
	std::vector<uint64_t> numNodesInLayer;
	uint64_t nodeCounter = 0;
	int graphDepth = -1;
	for (id_t nodeIndex : nodeOrder) {
		const auto &node = nodes.find(nodeIndex)->second;
		const int nodeDepth = node.depth;

		if (nodeDepth > graphDepth) {
			++graphDepth;
			m_nodes.push_back(std::vector<DrawableNode>());

			numNodesInLayer.push_back(nodeCounter);
			nodeCounter = 0;
		}

		++nodeCounter;
	}

	// Find scale
	// node spacing
	uint64_t maxNumNodesInAnyLayer = 0;
	for (uint64_t numNodes : numNodesInLayer) {
		if (numNodes > maxNumNodesInAnyLayer) {
			maxNumNodesInAnyLayer = numNodes;
		}
	}

	if (maxNumNodesInAnyLayer == 0) {
		m_nodeSpacing = 0;
	} else if (maxNumNodesInAnyLayer & 1) { // odd
		m_nodeSpacing = height / (maxNumNodesInAnyLayer - 1);
	} else {
		m_nodeSpacing = height / maxNumNodesInAnyLayer;
	}
	m_nodeSpacing *= 0.8;

	// layer spacing
	if (graphDepth == 0) {
		m_layerSpacing = 0;
	} else {
		m_layerSpacing = width * 0.7 / graphDepth;
	}

	// node radius
	m_nodeRadius = std::min(m_nodeSpacing / 5.0, m_layerSpacing / 5.0);

	// Add drawable nodes
	std::unordered_map<int, std::pair<int, int>> nodeIndexToPlace;
	for (id_t nodeIndex : nodeOrder) {
		const auto &node = nodes.find(nodeIndex)->second;
		const int nodeDepth = node.depth;

		addDrawableNode(nodeDepth, nodeIndex);

		const id_t nodeIndexAtDepth = m_nodes[nodeDepth].size() - 1;
		nodeIndexToPlace[nodeIndex] = std::make_pair(nodeDepth, nodeIndexAtDepth);
	}

//	int numConnections = 0;
	for (id_t nodeIndex : nodeOrder) {
		const auto &node = nodes.find(nodeIndex)->second;
		const auto &outputNodeIndices = node.outputNodeIds;

		const auto& nodePlace = nodeIndexToPlace[nodeIndex];
		const int depth = nodePlace.first;
		const int index = nodePlace.second;
		const sf::Vector2f& startPosition = m_nodes[depth][index].circle.getPosition();

		for (id_t outputNodeIndex : outputNodeIndices) {
			const auto& outputNodePlace = nodeIndexToPlace[outputNodeIndex];
			const int outputDepth = outputNodePlace.first;
			const int outputIndex = outputNodePlace.second;
			const sf::Vector2f& endPosition = m_nodes[outputDepth][outputIndex].circle.getPosition();

			m_connections.push_back(DrawableConnection(startPosition, endPosition, m_nodeRadius));

//			++numConnections;
		}
	}

//	std::cout << "Has loaded: " << numNodes << " nodes into: " << m_nodes.size() << " layers, with: " << numConnections << " connections!" << std::endl;
}

void DAGRenderer::render(sf::RenderWindow &window) {
	for (const auto& connection : m_connections) {
		window.draw(connection.arrowBody);
		window.draw(connection.arrowHead);
	}

	for (auto &nodesInLayer : m_nodes) {
		for (auto &node : nodesInLayer) {
			window.draw(node.circle);

			node.text.setFont(node.font);
			window.draw(node.text);
		}
	}
}
