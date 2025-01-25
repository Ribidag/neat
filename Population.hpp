#ifndef NEAT_POPULATION_HPP_
#define NEAT_POPULATION_HPP_

#include <SFML/Graphics.hpp>

#include "GenePool.hpp"
#include "rendering/ScoringVisualizer.hpp"
#include "rendering/RecordViewer.hpp"

class Phenotype;

class Population {
private:
	GenePool m_genePool;

	ScoringVisualizer m_scoringVisualizer;

	RecordViewer m_speciesRecordViewer;
	RecordViewer m_fitnessRecordViewer;

    uint64_t m_generationId;

    double m_topGenerationFitness;
    std::vector<double> m_fitnessRecord;
    id_t m_bestGenotypeIdInGeneration;

public:
    static void scorePhenotypeShell(std::unordered_map<id_t, std::unique_ptr<Phenotype>>::iterator phenotypeIt, const std::unordered_map<id_t, size_t>& vectorIndices, std::vector<double>& genotypeScoresVector, const uint32_t iterationsThisThread);

public:
	Population(const uint64_t numGenotypes, const uint16_t numInputs, const uint16_t numOutputs);

	static double scorePhenotype(const std::unique_ptr<Phenotype>& phenotypePtr);
	void scoreGenerationThreaded();
	void scoreGeneration();
	void selection();

	void perform(sf::RenderWindow& window, bool& showFitnessRecord, bool& showSpeciesRecord, bool& showGenotype, bool& showPerformance, double& time_step);

};

#endif /* NEAT_POPULATION_HPP_ */
