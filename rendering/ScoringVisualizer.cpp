#include "ScoringVisualizer.hpp"

#include "mvc/Model.hpp"
#include "mvc/View.hpp"

#include "../Genotype.hpp"

#include <iostream>

ScoringVisualizer::ScoringVisualizer(const std::array<double, 4>& boundingRect)
	: m_font()
	, m_width(0)
	, m_height(1)
	, m_dagRenderer()
	, m_scoreboardView()
	, m_genotypeView()
	, m_simulationView() {

	const double x0 = boundingRect[0];
	const double y0 = boundingRect[1];
	m_width = boundingRect[2];
	m_height = boundingRect[3];

	if (!m_font.loadFromFile("C:\\Windows\\Fonts\\arial.ttf")) {
		std::cerr << "NO FONT" << std::endl;
	}

	m_scoreboardView.setViewport(sf::FloatRect(x0 + 0.3 * m_width, y0, 0.7 * m_width, 0.2 * m_height));
	m_genotypeView.setViewport(sf::FloatRect(x0, y0, 0.3 * m_width, m_height));
	m_simulationView.setViewport(sf::FloatRect(x0 + 0.3 * m_width, y0 + 0.2*m_height, 0.7 * m_width, 0.8 * m_height));

	m_simulationView.setSize(1, 1);
}

void ScoringVisualizer::renderScoreboard(sf::RenderWindow &window, const double score) {
	const double windowWidth = window.getSize().x;
    const double windowHeight = window.getSize().y;
    m_scoreboardView.reset(sf::FloatRect(0, 0, windowWidth * 0.7 * m_width, windowHeight * 0.2 * m_height));

    const double width = m_scoreboardView.getSize().x;
    const double height = m_scoreboardView.getSize().y;
	sf::Text scoreText;
	scoreText.setString("Score " + std::to_string(score));
	scoreText.setCharacterSize(static_cast<int>(std::min(height / 2.0, width / 10.0)));
	scoreText.setFont(m_font);
	sf::Vector2f position(width / 10.0, height / 10.0);
	scoreText.setPosition(position);
	scoreText.setFillColor(sf::Color(200, 200, 200));

	window.setView(m_scoreboardView);
	window.draw(scoreText);
}

void ScoringVisualizer::loadGenotype(sf::RenderWindow &window, const Genotype &genotype) {
	const double windowWidth = window.getSize().x;
    const double windowHeight = window.getSize().y;
    m_genotypeView.reset(sf::FloatRect(0, 0, windowWidth * 0.3 * m_width, windowHeight * m_height));

    const float width = m_genotypeView.getSize().x;
    const float height = m_genotypeView.getSize().y;
	m_dagRenderer.loadDAG(genotype, width, height);
}

void ScoringVisualizer::renderGenotype(sf::RenderWindow &window) {
	window.setView(m_genotypeView);
	m_dagRenderer.render(window);
}

void ScoringVisualizer::renderSimulation(sf::RenderWindow &window, const Model &model) {
	const double windowWidth = window.getSize().x;
    const double windowHeight = window.getSize().y;
	m_simulationView.reset(sf::FloatRect(0, 0, windowWidth * 0.7 * m_width, windowHeight * 0.8 * m_height));
	View::render(model, window, m_simulationView);
}

//void ScoringVisualizer::renderSimulation(sf::RenderWindow &window, const FBModel &model) {
//	const double windowWidth = window.getSize().x;
//    const double windowHeight = window.getSize().y;
//	m_simulationView.reset(sf::FloatRect(0, 0, windowWidth * 0.7 * m_width, windowHeight * 0.8 * m_height));
//	FBView::render(model, window, m_simulationView);
//}
