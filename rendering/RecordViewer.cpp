#include "RecordViewer.hpp"

RecordViewer::RecordViewer(float x0, float y0, float x1, float y1)
	: m_view() {
	m_view.setViewport(sf::FloatRect(x0, y0, x1, y1));
}

void RecordViewer::renderSpeciesRecord(sf::RenderWindow& window, const std::vector<std::vector<double>>& speciesRecord) const {
    // Calculate the total height of each row
	const float rowHeight = m_view.getSize().y;

    // Calculate row widths
    const float viewWidth = m_view.getSize().x;
    const float rowWidth = viewWidth / speciesRecord.size();

    // Create a vertex array for drawing
    sf::VertexArray bars(sf::Quads);

    // Set up colors
    sf::Color color1 = sf::Color::Blue;
    sf::Color color2 = sf::Color::Red;

    // Position and draw the bars
    float xPosition = 0.0f;
    float yPosition = 0.0f;
    for (const auto& row : speciesRecord) {
        float yOffset = 0.0f;
        for (size_t i = 0; i < row.size(); ++i) {
            float height = row[i] * rowHeight; // Scale the height
            sf::Color color = (i % 2 == 0) ? color1 : color2;

            // Add vertices for the bar
            bars.append(sf::Vertex(sf::Vector2f(xPosition, yPosition + yOffset), color));
            bars.append(sf::Vertex(sf::Vector2f(xPosition + rowWidth, yPosition + yOffset), color));
            bars.append(sf::Vertex(sf::Vector2f(xPosition + rowWidth, yPosition + yOffset + height), color));
            bars.append(sf::Vertex(sf::Vector2f(xPosition, yPosition + yOffset + height), color));

            yOffset += height;
        }
        xPosition += rowWidth; // Spacing between bars
    }

    window.setView(m_view);
    window.draw(bars);
}

//void RecordViewer::renderSpeciesRecord(sf::RenderWindow& window, const std::vector<std::vector<double>>& speciesRecord) const {
//    // Calculate the total height of each row
//	const float rowHeight = m_view.getSize().y;
//
//    std::vector<double> rowHeights;
//    for (const auto& row : speciesRecord) {
//        double totalHeight = 0.0;
//        for (double portion : row) {
//            totalHeight += portion;
//        }
//        rowHeights.push_back(totalHeight);
//    }
//
//    // Calculate row widths
//    const float viewWidth = m_view.getSize().x;
//    const float rowWidth = viewWidth / speciesRecord.size();
//
//    // Create a vertex array for drawing
//    sf::VertexArray bars(sf::Quads);
//
//    // Set up colors
//    sf::Color color1 = sf::Color::Blue;
//    sf::Color color2 = sf::Color::Red;
//
//    // Position and draw the bars
//    float xPosition = 0.0f;
//    float yPosition = 0.0f;
//    for (const auto& row : speciesRecord) {
//        float yOffset = 0.0f;
//        for (size_t i = 0; i < row.size(); ++i) {
//            float height = row[i] * rowHeights[i] * rowHeight; // Scale the height
//            sf::Color color = (i % 2 == 0) ? color1 : color2;
//
//            // Add vertices for the bar
//            bars.append(sf::Vertex(sf::Vector2f(xPosition, yPosition + yOffset), color));
//            bars.append(sf::Vertex(sf::Vector2f(xPosition + rowWidth, yPosition + yOffset), color));
//            bars.append(sf::Vertex(sf::Vector2f(xPosition + rowWidth, yPosition + yOffset + height), color));
//            bars.append(sf::Vertex(sf::Vector2f(xPosition, yPosition + yOffset + height), color));
//
//            yOffset += height;
//        }
//        xPosition += rowWidth; // Spacing between bars
//    }
//
//    window.setView(m_view);
//    window.draw(bars);
//}

void RecordViewer::renderFitnessRecord(sf::RenderWindow &window, const std::vector<double>& fitnessRecord) const {
    sf::VertexArray line(sf::LineStrip);

    const float viewHeight = m_view.getSize().y;
    const float middleHeight = viewHeight / 2.0f;
    const float viewWidth = m_view.getSize().x;
    const float width = viewWidth / fitnessRecord.size();

    float max = *std::max_element(fitnessRecord.begin(), fitnessRecord.end());

    for (size_t i = 0; i < fitnessRecord.size(); ++i) {
        const float x = static_cast<float>(i * width);

        float y = static_cast<float>(fitnessRecord[i]);
        y = y/max * viewHeight;

        y = 2*middleHeight - y; // flip

        line.append(sf::Vertex(sf::Vector2f(x, y), sf::Color::Blue));
    }

    window.setView(m_view);
    window.draw(line);
}
