#ifndef NEAT_RENDERING_SCORINGVISUALIZER_HPP_
#define NEAT_RENDERING_SCORINGVISUALIZER_HPP_

#include <SFML/Graphics.hpp>

#include "DAGRenderer.hpp"

class Model;
class FBModel;
class Genotype;

class ScoringVisualizer {
private:
	sf::Font m_font;

	double m_width;
	double m_height;

	DAGRenderer m_dagRenderer;

	sf::View m_scoreboardView;
	sf::View m_genotypeView;
	sf::View m_simulationView;

public:
	ScoringVisualizer(const std::array<double, 4>& boundingRect);

	void loadGenotype(sf::RenderWindow &window, const Genotype& genotype);

	void renderScoreboard(sf::RenderWindow& window, const double score);
	void renderGenotype(sf::RenderWindow& window);
	void renderSimulation(sf::RenderWindow& window, const Model& model);
//	void renderSimulation(sf::RenderWindow& window, const FBModel& model);

	void renderScoring(sf::RenderWindow& window, const Genotype& genotype);

};

#endif /* NEAT_RENDERING_SCORINGVISUALIZER_HPP_ */
