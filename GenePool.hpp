#ifndef NEAT_GENEPOOL_HPP_
#define NEAT_GENEPOOL_HPP_

#include <unordered_map>
#include <map>

#include "util/types.hpp"
#include "Genotype.hpp"

#include <SFML/Graphics.hpp>

class GenePool {
private:
	struct Species {
		std::vector<id_t> genotypeIds;
		Genotype representativeGenotype;

		id_t championId;
		double adjustedFitness;

		double speciesGenotypeFitnessRecord = 0;
		uint64_t generationsWithoutImprovement = 0;
	};

	static constexpr double c1 = 2;
	static constexpr double c2 = 2;
	static constexpr double c3 = 3.0; // 0.4 | 3.0

	static constexpr double geneticDistanceBoundary = 4.0; // 3.0 | 4.0
	static constexpr double interspeciesCrossoverProbability = 0.001;
	static constexpr double noCrossoverProbability = 0.25;
	static constexpr double splitSynapseProbability = 0.03; // 0.03
	static constexpr double growSynapseProbability = 0.3; // 0.05 | 0.3
	static constexpr double mutateSynapseWeightProbability = 0.8;
	static constexpr double shiftSynapseWeightProbability = 0.9;
	static constexpr double mutateNeuronBiasProbability = 0.5;
	static constexpr double shiftNeuronBiasProbability = 0.95;

	const uint64_t m_numGenotypes;
	const uint16_t m_numInputs;
	const uint16_t m_numOutputs;

	innov_t m_innovationIndex;
	innov_t m_neuronIndex;

	id_t m_speciesIndex;
	std::map<id_t, Species> m_species;
	std::vector<id_t> m_speciesIds;
	std::unordered_map<id_t, uint64_t> m_speciesNumberOfNextGeneration;
	std::unordered_set<id_t> m_championIds;
	std::vector<std::vector<double>> m_speciesRecord;

	id_t m_genotypeIndex;
	std::unordered_map<id_t, Genotype> m_genotypes;
	std::unordered_map<id_t, double> m_genotypeScores;

	uint64_t m_generationId;

	double m_genotypeFitnessRecord;

public:
	static double findGeneticDistance(const Genotype& genotype1, const Genotype& genotype2);

private:
	const id_t getRandomSpeciesId() const;
	const std::vector<double> getCumulativeFitness(const std::vector<id_t>& genotypeIds) const;
	const id_t getRandomGenotypeId(const std::vector<id_t>& genotypeIds, const std::vector<double> &cumulativeFitness) const;

	void speciate();
	void removeWeakerGenotypes();
	void removeWeakerSpecies();
	void assignOffspringToSpecies();
	void mateGenotypes();
	void mutateGenotypes();

public:
	GenePool(const uint64_t numGenotypes, const uint16_t numInputs, const uint16_t numOutputs);

	void constructNextGeneration();

	const std::unordered_map<id_t, Genotype>& getGenotypes() const {
		return m_genotypes;
	}

	std::unordered_map<id_t, Genotype>& getGenotypes() {
		return m_genotypes;
	}

	void setGenotypeScores(const std::unordered_map<id_t, double> &genotypeScores) {
		m_genotypeScores = genotypeScores;
	}

	double getGenotypeFitnessRecord() const {
		return m_genotypeFitnessRecord;
	}

public:
	const std::vector<std::vector<double> >& getSpeciesRecord() const {
		return m_speciesRecord;
	}
};

#endif /* NEAT_GENEPOOL_HPP_ */
