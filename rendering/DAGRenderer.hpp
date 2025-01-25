#ifndef NEAT_RENDERING_DAGRENDERER_HPP_
#define NEAT_RENDERING_DAGRENDERER_HPP_

#include <vector>
#include <cmath>

#include <SFML/Graphics.hpp>

#include "../util/types.hpp"

#include <iostream>

class DirectedAcyclicGraph;

class DAGRenderer {
private:
	struct DrawableNode {
		sf::CircleShape circle;
		sf::Font font;
		sf::Text text;

		DrawableNode(const sf::Vector2f& position, const double radius, const id_t id)
			: circle(radius, 20), font(), text() {

			circle.setOrigin(radius, radius);
			circle.setPosition(position);

			if (!font.loadFromFile("C:\\Windows\\Fonts\\arial.ttf")) {
				std::cerr << "NO FONT" << std::endl;
			}
			text.setString(std::to_string(id));
			text.setCharacterSize(static_cast<int>(radius * 1.5));
			text.setFont(font);
			sf::FloatRect textRect = text.getLocalBounds();
			text.setOrigin(textRect.left + textRect.width/2.0f, textRect.top  + textRect.height/2.0f);
			text.setPosition(position);
			text.setFillColor(sf::Color(0, 0, 0));
		}
	};

	struct DrawableConnection {
		sf::RectangleShape arrowBody;
		sf::CircleShape arrowHead;

		DrawableConnection(const sf::Vector2f &startPosition, const sf::Vector2f &endPosition, const float nodeRadius) {
			float x0, y0, x1, y1, dx, dy, width, height;
			x0 = startPosition.x;
			y0 = startPosition.y;
			x1 = endPosition.x;
			y1 = endPosition.y;

			dx = x1 - x0;
			dy = y1 - y0;

			double angle = std::atan2(dy, dx)  * 180.0 / M_PI;

			width = std::sqrt(dx*dx + dy*dy);
			height = 1 * nodeRadius / 5.0;

			arrowBody.setSize(sf::Vector2f(width, height));
			arrowBody.setOrigin(0, 0.5*height);
			arrowBody.setPosition(x0, y0);
			arrowBody.setRotation(angle);

			arrowHead.setPointCount(3);
			const float radius = 2 * nodeRadius / 5.0;
			arrowHead.setRadius(radius);
			arrowHead.setOrigin(radius, radius);
			arrowHead.setRotation(angle + 90);
			arrowHead.setPosition(startPosition + (endPosition - startPosition) * (width - nodeRadius - radius*0.6f) / width);
		}
	};

	sf::Vector2f m_position;

	float m_nodeSpacing;
	float m_layerSpacing;
	float m_nodeRadius;

	std::vector<std::vector<DrawableNode>> m_nodes;
	std::vector<DrawableConnection> m_connections;

private:
	void addDrawableNode(const uint32_t depth, const id_t id);

public:
	DAGRenderer();

	void loadDAG(const DirectedAcyclicGraph &dag, const float width, const float height);
	void render(sf::RenderWindow& window);
};

#endif /* NEAT_RENDERING_DAGRENDERER_HPP_ */
