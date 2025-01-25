#include <iostream>
#include <SFML/Graphics.hpp>
#include "util/NumberGenerator.hpp"
#include "Interface.hpp"

#include "Population.hpp"

int main() {
	NumberGenerator::initialize();
//	Interface::initialize();

	const int width = 800;
	const int height = 600;
	sf::RenderWindow window(sf::VideoMode(width, height), "NEAT");

	Population population(200, 7, 4); // Snake : (8, 3), Flappy bird : (3, 1) // Drone : (7, 4)

	for (int i(0); i < 0; ++i) {
		NumberGenerator::getSym(1.0);
	}

	bool showFitnessRecord = true;
	bool showSpeciesRecord = false;
	bool showGenotype = true;
	bool showPerformance = true;
	double msPerFrame = 10;

	while (window.isOpen()) {
		sf::Event event;
		while (window.pollEvent(event)) {
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				if (event.key.code == sf::Keyboard::P) {
					showPerformance = !showPerformance;
				} else if (event.key.code == sf::Keyboard::Escape)
					window.close();
				break;
			default:
				break;
			}
		}

//		std::cout << "Scoring Generation!!" << std::endl;
		population.scoreGenerationThreaded();
//		std::cout << "Scored Generation!!" << std::endl;
		if (showPerformance) {
			population.perform(window, showFitnessRecord, showSpeciesRecord, showGenotype, showPerformance, msPerFrame);
		}
		population.selection();
	}

	return 0;
}
