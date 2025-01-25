#ifndef NEAT_RENDERING_RECORDVIEWER_HPP_
#define NEAT_RENDERING_RECORDVIEWER_HPP_

#include <vector>

#include <SFML/Graphics.hpp>

#include "../util/types.hpp"

class RecordViewer {
private:
	sf::View m_view;

public:
	RecordViewer(float x0, float y0, float x1, float y1);

	void renderSpeciesRecord(sf::RenderWindow& window, const std::vector<std::vector<double>>& speciesRecord) const;
	void renderFitnessRecord(sf::RenderWindow& window, const std::vector<double>& fitnessRecord) const;
};

#endif /* NEAT_RENDERING_RECORDVIEWER_HPP_ */
