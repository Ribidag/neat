#include "GenePool.hpp"

#include <cmath>
#include <algorithm>

#include "util/NumberGenerator.hpp"

#include <iostream>
#include "rendering/DAGRenderer.hpp"

GenePool::GenePool(const uint64_t numGenotypes, const uint16_t numInputs, const uint16_t numOutputs)
	: m_numGenotypes(numGenotypes)
	, m_numInputs(numInputs)
	, m_numOutputs(numOutputs)
	, m_innovationIndex(numInputs * numOutputs)
	, m_neuronIndex(numOutputs)
	, m_speciesIndex(0)
	, m_species()
	, m_speciesIds()
	, m_speciesNumberOfNextGeneration()
	, m_championIds()
	, m_speciesRecord()
	, m_genotypeIndex(numGenotypes-1)
	, m_generationId(0)
	, m_genotypeFitnessRecord(0) {

	for (uint64_t i(0); i < numGenotypes; ++i) {
		m_genotypes[i] = Genotype(numInputs, numOutputs);
	}
}

double GenePool::findGeneticDistance(const Genotype &genotype1, const Genotype &genotype2) {

	const auto& synapseGenes1 = genotype1.getSynapseGenes();
	const auto& synapseGenes2 = genotype2.getSynapseGenes();
	const innov_t latestInnovation1 = genotype1.getLatestInnovation();
	const innov_t latestInnovation2 = genotype2.getLatestInnovation();

	innov_t maxInnovation, edgeInnovation;
	std::unordered_map<innov_t, SynapseGene> derivedGenes, basalGenes;
	if (latestInnovation1 > latestInnovation2) {
		maxInnovation = latestInnovation1;
		edgeInnovation = latestInnovation2;

		derivedGenes = synapseGenes1;
		basalGenes = synapseGenes2;
	} else {
		maxInnovation = latestInnovation2;
		edgeInnovation = latestInnovation1;

		derivedGenes = synapseGenes2;
		basalGenes = synapseGenes1;
	}

	uint32_t numMatchingGenes = 0;
	uint32_t numDisjointGenes = 0;
	uint32_t numExcessGenes = 0;
	double weightDifferenceAverage;
	uint32_t normalizationConstant;

	if (synapseGenes1.size() >= synapseGenes2.size()) {
		normalizationConstant = synapseGenes1.size();
	} else {
		normalizationConstant = synapseGenes2.size();
	}
	if (normalizationConstant < 20) {
		normalizationConstant = 1;
	}

	double weightDifferenceSum = 0;
	for (innov_t i(1); i <= maxInnovation; ++i) {
		const auto& derivedIt = derivedGenes.find(i);

		if (i <= edgeInnovation) {
			const auto& basalIt = basalGenes.find(i);

			// matching
			if (derivedIt != derivedGenes.end() && basalIt != basalGenes.end()) {
				++numMatchingGenes;

				weightDifferenceSum += std::abs(derivedIt->second.getWeight() - basalIt->second.getWeight());

			} else if (derivedIt != derivedGenes.end() || basalIt != basalGenes.end()) {
				++numDisjointGenes;
			}

		} else if (derivedIt != derivedGenes.end()) {
				++numExcessGenes;
		}
	}

	weightDifferenceAverage = weightDifferenceSum / numMatchingGenes;

	//std::cout << "		Excess: " << numExcessGenes << ", Disjoint: " << numDisjointGenes << ", W: " << weightDifferenceAverage << ", Matching: " << numMatchingGenes << std::endl;

	double biasDifferenceSum = 0;
	for (const auto& neuronGeneIt1 : genotype1.getNeuronGenes()) {
		const id_t neuronGeneId = neuronGeneIt1.first;
		const double bias1 = neuronGeneIt1.second;

		const auto& neuronGeneIt2 = genotype2.getNeuronGenes().find(neuronGeneId);
		if (neuronGeneIt2 != genotype2.getNeuronGenes().end()) {
			const double bias2 = neuronGeneIt2->second;
			biasDifferenceSum += std::abs(bias1 - bias2);
		}
	}

	double biasDifferenceAverage = biasDifferenceSum / numMatchingGenes;

	return (c1 * numExcessGenes + c2 * numDisjointGenes) / normalizationConstant + c3 * (weightDifferenceAverage + biasDifferenceAverage);
}

const id_t GenePool::getRandomSpeciesId() const {
	uint64_t total = 0;
	const uint64_t randomRoll = std::floor(NumberGenerator::getASym(m_numGenotypes));

	for (const id_t id : m_speciesIds) {
		const auto& species = m_species.find(id)->second;
		uint64_t numGenotypesInSpecies = species.genotypeIds.size();

		if (randomRoll >= total && randomRoll < total + numGenotypesInSpecies) {
			return id;
		}

		total += numGenotypesInSpecies;
	}

	std::cerr << "Failed to find random speciesId" << std::endl;
	return 0;
}

const std::vector<double> GenePool::getCumulativeFitness(const std::vector<id_t> &genotypeIds) const {
	double minimumFitness = -1;
	for (const id_t genotypeId : genotypeIds) {
		const double genotypeFitness = m_genotypeScores.find(genotypeId)->second;

		if (genotypeFitness < minimumFitness || minimumFitness < 0) {
			minimumFitness = genotypeFitness;
		}
	}

	double totalFitness = 0.0;
	for (const id_t genotypeId : genotypeIds) {
		const double genotypeFitness = m_genotypeScores.find(genotypeId)->second - minimumFitness;
		totalFitness += genotypeFitness;
	}

//	std::cout << "Shares: ";
	std::vector<double> cumulativeFitness;
	double partialFitness = 0.0;
	for (const id_t genotypeId : genotypeIds) {
		const double genotypeFitness = m_genotypeScores.find(genotypeId)->second - minimumFitness;
		const double genotypeShare = genotypeFitness / totalFitness;
		partialFitness += genotypeShare;
		cumulativeFitness.push_back(partialFitness);
//		std::cout << genotypeShare << " (" << genotypeFitness << ") ";
	}
//	std::cout << std::endl;

	return cumulativeFitness;
}


const id_t GenePool::getRandomGenotypeId(const std::vector<id_t> &genotypeIds, const std::vector<double> &cumulativeFitness) const {
	double randomValue = NumberGenerator::getASym(1.0);

	for (size_t i = 0; i < cumulativeFitness.size(); ++i) {
		if (randomValue < cumulativeFitness[i]) {
			return genotypeIds[i];
		}
	}

	// Fallback: return the first element
	return genotypeIds[0];
}

void GenePool::speciate() {
//	std::cout << "Speciating!" << std::endl;

	/*
	 * For every genotype, go through each species and check the genetic distance between
	 * the genotype and that species' representative genotype. If the genetic distance is
	 * within the genetic distance boundary, assign the genotype to that species and stop
	 * the search.
	 *
	 * If no matching species is found, a new species is created with that genotype as its
	 * representative genotype.
	 *
	 * Some of the old species might not recieve any genotypes at all. In that case, they
	 * will be deleted later, as they have gone extinct.
	 */
	for (const auto& genotypeIt : m_genotypes) {
		const id_t genotypeId = genotypeIt.first;
		const auto& genotype = genotypeIt.second;

//		std::cout << "Attempting to determine species of genotype " << genotypeId << std::endl;

		/*
		 * Here is the loop to go through all species.
		 */
		bool speciesFound = false;
		for (auto& speciesIt : m_species) {
			// Find representativeGenotype of the species.
			auto& species = speciesIt.second;
			const auto& representativeGenotype = species.representativeGenotype;

//			auto& speciesId = speciesIt.first;
//			std::cout << "	Comparing to species " << speciesId << std::endl;

			// Assign the genotype to the species if the geneticDistance is within the threshold.
			const double geneticDistance = findGeneticDistance(genotype, representativeGenotype);
			if (geneticDistance <= geneticDistanceBoundary) {
				speciesFound = true;
				species.genotypeIds.push_back(genotypeId);

//				std::cout << "		Match!" << std::endl;
				break;
			}
		}

		/*
		 * If no matching species is found, create a new one with the genotype
		 * as its representative.
		 *
		 * Note here that the species index is assigned and then incremented.
		 */
		if (!speciesFound) {
//			std::cout << "	No species found! Creating species " << m_speciesIndex << std::endl;

			Species species;
			species.genotypeIds.push_back(genotypeId);
			species.representativeGenotype = genotype;
			m_species[m_speciesIndex] = species;
			m_speciesIds.push_back(m_speciesIndex);
			++m_speciesIndex;
		}
	}

//	for (auto &speciesIt : m_species) {
//		const id_t speciesId = speciesIt.first;
//		auto &species = speciesIt.second;
//		auto &genotypeIds = species.genotypeIds;
//		const uint64_t numGenotypesInSpecies = genotypeIds.size();
//
//		std::cout << "Species " << speciesId << " has " << numGenotypesInSpecies << " genotypes:" << std::endl;
//
//		for (const id_t genotypeId : genotypeIds) {
//			const auto& genotype = m_genotypes[genotypeId];
//			std::cout << "	" << genotypeId << ": Genetic distance " << findGeneticDistance(genotype, species.representativeGenotype) << " to species representative" << std::endl;
//		}
//	}
}

void GenePool::removeWeakerGenotypes() {
	/*
	 * Each species will firstly sum the fitnesses of its genotypes to determine the
	 * species fitness. Then, this total species fitness is divided by the
	 * number of genotypes in the species to obtain the adjusted species fitness.
	 * The adjusted species fitness is what is used to determine how big a share of
	 * the next generation this species is entitled to. Because larger species are
	 * penalized, there is room for less dominant species to explore other niches.
	 *
	 * Additionally, the highest fitness from any genotype in the current generation
	 * of the species is noted.
	 *
	 * Then, the species will order its genotypes from best to worst performing
	 * and delete the bottom half of genotypes (there will always be at least one
	 * genotype left). Thus, genotypes compete only with other genotypes within the
	 * same species.
	 *
	 * The species will also designate its champion. If a species is large enough, its
	 * champion will later be copied directly into the next generation.
	 *
	 * Lastly, the species compares the top genotype fitness in this generation to its
	 * record over all past generations. If the record is not beat, the species counts
	 * one more generation without improvement. Species which do not improve over many
	 * generations will later be eliminated.
	 */
	for (auto &speciesIt : m_species) {
		// Species info
		auto &species = speciesIt.second;
		auto &genotypeIds = species.genotypeIds;
		const uint64_t numGenotypesInSpecies = genotypeIds.size();

		/*
		 * Here the adjusted species fitness and the top genotype fitness are calculated.
		 */
		double speciesFitness = 0;
		double topFitnessInGeneration = 0;

		for (const id_t genotypeId : genotypeIds) {
			const double genotypeFitness = m_genotypeScores[genotypeId];
			speciesFitness += genotypeFitness;

			topFitnessInGeneration = std::max(topFitnessInGeneration, genotypeFitness);
			m_genotypeFitnessRecord = std::max(m_genotypeFitnessRecord, genotypeFitness);
		}

		species.adjustedFitness = speciesFitness / numGenotypesInSpecies;

		//	const id_t speciesId = speciesIt.first;
		//	std::cout << "Species " << speciesId << " has unadjusted fitness: " << unadjustedSpeciesFitness << " from " << numGenotypesInSpecies << " genotypes." << std::endl;

		/*
		 * Here we sort the genotypes from highest fitness to lowest fitness.
		 * The best genotype is declared the champion of the species. Then,
		 * the worst performing 50% are eliminated.
		 */
		std::sort(genotypeIds.begin(), genotypeIds.end(), [this](const id_t i, const id_t j) {
		    return m_genotypeScores[i] > m_genotypeScores[j];
		});

		// Crown species champion as the first genotype in the list.
		species.championId = genotypeIds[0];

//		std::cout << "Scored order: ";
//		for (const id_t id : genotypeIds) {
//			std::cout << id << " (" << m_genotypeScores[id] << ")" << "	";
//		}
//		std::cout << std::endl;

		const size_t keepingIndex = std::ceil(numGenotypesInSpecies / 2.0);
		const int genotypesToRemove = numGenotypesInSpecies - keepingIndex;

//		std::cout << "	Erasing genotypes: ";
		for (int i(0); i < genotypesToRemove; ++i) {
			const id_t genotypeId = genotypeIds.back();
			genotypeIds.pop_back();
			m_genotypes.erase(genotypeId);
			m_genotypeScores.erase(genotypeId);
//			std::cout << genotypeId << " ";
		}
//		std::cout << std::endl;

		/*
		 * Here, the species tracks its current fitness in comparison to
		 * its best performance.
		 */
		if (topFitnessInGeneration <= species.speciesGenotypeFitnessRecord) {
			++species.generationsWithoutImprovement;
		} else {
			species.generationsWithoutImprovement = 0;
			species.speciesGenotypeFitnessRecord = topFitnessInGeneration;
		}
	}

	/*
	 * Species record.
	 */
	m_speciesRecord.emplace_back();
	m_speciesRecord.back().reserve(m_species.size());
	for (auto &speciesIt : m_species) {
		auto &species = speciesIt.second;
		const uint64_t numGenotypesInSpecies = species.genotypeIds.size();

		m_speciesRecord.back().push_back(static_cast<double>(numGenotypesInSpecies) / static_cast<double>(m_genotypes.size()));
	}
}

template <typename T>
int signum(T val) {
    return (T(0) < val) - (val < T(0));
}

void GenePool::removeWeakerSpecies() {
//	std::cout << "Erasing species: ";

	/*
	 * Those species with a single genotype and those which have not improved
	 * for 15 generations will be deleted. This is to save space for those species
	 * which seem to be progressing. There is now the possibility that all species
	 * will die out, resulting in the failure of the program.
	 */
	const size_t startingNumOfSpecies = m_speciesIds.size();
	size_t removed = 0;
	for (size_t i(0); i < startingNumOfSpecies - removed;) {
		const id_t speciesId = m_speciesIds[i];

		const bool speciesIsExtinct = m_species[speciesId].genotypeIds.size() == 0;
		const bool speciesIsNotImproving = m_species[speciesId].generationsWithoutImprovement >= 15;
		bool speciesIsNotImprovingAndBad = false;

		if (speciesIsNotImproving) {
			const double topGenerationFitness = m_species[speciesId].speciesGenotypeFitnessRecord;

			if (topGenerationFitness < 0.9*m_genotypeFitnessRecord) {
				speciesIsNotImprovingAndBad = true;
			}
		}

		// || m_species[speciesId].generationsWithoutImprovement >= 15
		if (speciesIsExtinct || speciesIsNotImprovingAndBad) {

//			std::cout << "Deleted a species with top fitness: " << m_species[speciesId].speciesGenotypeFitnessRecord << " / " << m_genotypeFitnessRecord << " and " << m_species[speciesId].genotypeIds.size() << " genotypes." << std::endl;

			std::swap(m_speciesIds[i], m_speciesIds.back());
			m_speciesIds.pop_back();
			m_species.erase(speciesId);
			++removed;

//			std::cout << speciesId << " ";
		} else {
			++i;
		}
	}
//	std::cout << std::endl;

	if (m_species.size() == 0) {
		std::cerr << "Extinction!" << std::endl;
	}
}

void GenePool::assignOffspringToSpecies() {
	/*
	 * Each species will be awarded a share of the next generation in proportion
	 * to its contribution to the total adjusted fitness of the entire gene pool.
	 *
	 * However, because it is not possible to exactly split a discrete number of
	 * genotypes along percentage lines, the extra unassigned genotypes are allotted
	 * randomly.
	 */
	/*
	 * Here the total gene pool adjusted fitness is calculated.
	 */
//	std::cout << "Species fitnesses: ";
	double genePoolAdjustedFitness = 0;
	for (const auto &speciesIt : m_species) {
//		const id_t speciesId = speciesIt.first;
		const auto &species = speciesIt.second;
		genePoolAdjustedFitness += species.adjustedFitness;
//		std::cout << species.totalFitness << " (" << speciesId << ") ";
	}
//	std::cout << std::endl;

	/*
	 * Here a number of offspring are assigned to each species in accordance
	 * with its share of the total adjusted fitness.
	 */
//	std::cout << "Species offspring: ";
	m_speciesNumberOfNextGeneration.clear();
	for (const auto &speciesIt : m_species) {
		const id_t speciesId = speciesIt.first;
		const auto &species = speciesIt.second;

		const uint64_t number = std::floor(species.adjustedFitness / genePoolAdjustedFitness * m_numGenotypes);
		m_speciesNumberOfNextGeneration[speciesId] = number;
//		std::cout << number << " (" << speciesId << ") ";
	}
//	std::cout << std::endl;

	/*
	 * Here the last few left over genotypes are allotted randomly.
	 */
	uint64_t unadjustedNumber = 0;
	for (const auto &numberIt : m_speciesNumberOfNextGeneration) {
		const uint64_t number = numberIt.second;
		unadjustedNumber += number;
	}
	const int numberDifference = m_numGenotypes - unadjustedNumber;
//	std::cout << "Has crudely alloted " << unadjustedNumber << " genotypes to the next generation!" << std::endl;

	auto it = m_speciesNumberOfNextGeneration.begin();
	for (int i(0); i < std::abs(numberDifference); ++i) {
		++it->second;
		it++;
	}

	//TODO: fixa om ingen förbättring på 20 generationer så får bara två arter fortplantas

//	uint64_t adjustedNumber = 0;
//	for (const auto &numberIt : m_speciesNumberOfNextGeneration) {
//		const uint64_t number = numberIt.second;
//		adjustedNumber += number;
//
//		const id_t speciesId = numberIt.first;
//		std::cout << "	Species " << speciesId << " has been allotted " << number << " offspring" << std::endl;
//	}
//	std::cout << "Has adjusted and allotted " << adjustedNumber << " genotypes to the next generation!" << std::endl;
//
//	if (adjustedNumber != m_numGenotypes) {
//		std::cerr << "Population size not maintained!!" << std::endl;
//	}
}

void GenePool::mateGenotypes() {
	/*
	 * Every species fills its allotted number of offspring to construct the genotypes
	 * in the next generation.
	 *
	 * Before this, however, the species champions from species with 5 or more genotypyes
	 * are copied directly into the next generation. Their ids are also stored to spare
	 * them from later mutations.
	 *
	 * The remaining allotted slots for each species are then to be filled from the offspring
	 * of parent genotypes in this generation. The probability that a genotype is selected as
	 * a parent is proportional to its share of the total fitness of its species.
	 *
	 * A first parent is selected as described above. There are now three possibilities as
	 * for how an offspring is created.
	 * Firstly, if there is no crossover, the first parent genotype is copied directly into
	 * the next generation. Importantly, it can still be mutated later.
	 * Secondly, if there is interspecies crossover, a random species is selected and a second
	 * parent chosen with probability in proportion to its share of that other species' total
	 * fitness. The first and second parents are then crossed over to produce a genotype offspring
	 * for the first species.
	 * Thirdly, and most commonly, there is regular intraspecies crossover. Here the second parent
	 * is selected from the same species and crossed over with the first parent.
	 *
	 * Lastly, each species chooses a new species representative entirely at random. Note that
	 * these representatives therefore come from the current ("old") generation. The genotypes
	 * of the next generation will not be assigned to any species until after their evaluation.
	 */
	std::unordered_map<id_t, Genotype> nextGenotypes;

	/*
	 * Here we copy over each champion unchanged.
	 */
	m_championIds.clear();
	for (const auto &speciesIt : m_species) {
		const id_t speciesId = speciesIt.first;
		const Species& species = speciesIt.second;

		// species.genotypeIds.size() >= 5 &&
		if (species.genotypeIds.size() >= 5 && m_speciesNumberOfNextGeneration.find(speciesId)->second > 0) {
			const id_t speciesChampionId = species.championId;

			nextGenotypes[speciesChampionId] = m_genotypes[speciesChampionId];
			--m_speciesNumberOfNextGeneration.find(speciesId)->second;

			m_championIds.insert(speciesChampionId);
//			std::cout << "	Saving champion " << speciesChampionId << " with score: " << m_genotypeScores[speciesChampionId] << " from species " << speciesId << std::endl;
		}
	}

	/*
	 * Creating m_numGenotypes new genotypes through crossover
	 */
	for (const auto &numberIt : m_speciesNumberOfNextGeneration) {
		const uint64_t allottedGenotypesForSpecies = numberIt.second;

		const uint64_t firstParentSpeciesId = numberIt.first;
		const auto& firstParentSpecies = m_species[firstParentSpeciesId];

		const auto& firstParentSpeciesCumulativeFitness = getCumulativeFitness(firstParentSpecies.genotypeIds);

		for (uint64_t i(0); i < allottedGenotypesForSpecies; ++i) {
			++m_genotypeIndex;

			// Selecting the first parent genotype at random
			const uint64_t firstParentGenotypeId = getRandomGenotypeId(firstParentSpecies.genotypeIds, firstParentSpeciesCumulativeFitness);
			const auto& firstParentGenotype = m_genotypes[firstParentGenotypeId];
			const double firstParentFitness = m_genotypeScores[firstParentGenotypeId];

//			std::cout << "First parent is " << firstParentGenotypeId << " (" << m_genotypeScores[firstParentGenotypeId] << ") " <<" of species " << firstParentSpeciesId << ", ";

			// Selecting the second parent using a random crossover type
			const double rollForCrossover = NumberGenerator::getASym(1.0);

			// No crossover
			if (rollForCrossover >= 1.0 - noCrossoverProbability) {
//				std::cout << "no crossover";
				nextGenotypes[m_genotypeIndex] = firstParentGenotype;
			} else {

				id_t secondParentGenotypeId;

				// Interspecies crossover: choose any genotype for the second parent
				if (rollForCrossover < interspeciesCrossoverProbability) {
					const id_t secondParentSpeciesId = m_speciesIds[std::floor(NumberGenerator::getASym(m_speciesIds.size()))];
					const auto& secondParentSpecies = m_species[secondParentSpeciesId];

					const auto& secondParentSpeciesCumulativeFitness = getCumulativeFitness(secondParentSpecies.genotypeIds);
					secondParentGenotypeId = getRandomGenotypeId(secondParentSpecies.genotypeIds, secondParentSpeciesCumulativeFitness);

//					std::cout << "Second parent (interspecies) is " << secondParentGenotypeId << " of species " << secondParentSpeciesId;
				} // Intraspecies crossover: choose a genotype from the same species for the second parent
				else {
					secondParentGenotypeId = getRandomGenotypeId(firstParentSpecies.genotypeIds, firstParentSpeciesCumulativeFitness);
//					std::cout << "Second parent is " << secondParentGenotypeId << " of species " << firstParentSpeciesId;
				}

				const auto& secondParentGenotype = m_genotypes[secondParentGenotypeId];
				const double secondParentFitness = m_genotypeScores[secondParentGenotypeId];

				if (firstParentFitness > secondParentFitness) {
					nextGenotypes[m_genotypeIndex] = Genotype(firstParentGenotype, secondParentGenotype);
				} else {
					nextGenotypes[m_genotypeIndex] = Genotype(secondParentGenotype, firstParentGenotype);
				}
			}

//			std::cout << ". Roll was: " << rollForCrossover << std::endl;
//			std::cout << "	Created child " << m_genotypeIndex << std::endl;
		}
	}

	// set new species representatives
	for (auto& speciesIt : m_species) {
//		auto& speciesId = speciesIt.first;
		auto& species = speciesIt.second;
		if (species.genotypeIds.size() > 0) {
			const id_t representativeGenotypeId = species.genotypeIds[std::floor(NumberGenerator::getASym(species.genotypeIds.size()))];
			species.representativeGenotype = m_genotypes[representativeGenotypeId];

//			std::cout << "New representative for species " << speciesId << " is genotype " << representativeGenotypeId << std::endl;
			species.genotypeIds.clear();
		}
	}

	m_genotypes = nextGenotypes;
	m_genotypeScores.clear();
}

void GenePool::mutateGenotypes() {
	/*
	 * The genotypes of the new generation are mutated. There are two types of structural
	 * mutations that can be applied. An existing synapse can be split by introducing a new
	 * node in the middle, while new synapses can grow between existing nodes.
	 *
	 * When a new structural mutation is assigned to a genotype it is not immediately applied.
	 * This is because the same mutation can arise independently in several genotypes. For
	 * example, many genotypes can have the synapse with a certain innovation number. If two
	 * genotypes are mutated by splitting that synapse, we want the new node being created to
	 * get the same identifier in both genotypes.
	 *
	 * Thus, while determining mutations for the genotypes, those mutations that arise are stored
	 * along with a list of all genotypes that should have that mutation applied. Then, each
	 * unique mutation is applied to all associated genotypes with the same identifier for the
	 * new structure.
	 *
	 * After these structural mutations have been applied, every synapse in every genotype has
	 * a probability of shifting its weight a little, and every node its bias.
	 */

	std::unordered_map<innov_t,std::vector<id_t>> synapseSplitMutations;
	std::unordered_map<std::pair<id_t, id_t>, std::vector<id_t>, SzudzikHash> growSynapseMutations;

	// track all structural mutations
	for (auto& genotypeIt : m_genotypes) {
		const id_t genotypeId = genotypeIt.first;
		const auto& genotype = genotypeIt.second;

		// champions should be unchanged
		if (m_championIds.count(genotypeId)) {
//			std::cout << "Champion " << genotypeId << " is not mutated structurally" << std::endl;
			continue;
		}

		double rollForGrowSynapseMutation = NumberGenerator::getASym(1.0);
		double rollForSplitSynapseMutation = NumberGenerator::getASym(1.0);
//		std::cout << "Structural mutation step for new genotype " << genotypeId << " with rolls " << rollForGrowSynapseMutation << " (grow synapse) and " << rollForSplitSynapseMutation << " (split synapse)" << std::endl;

		// grow synapse mutation
		if (rollForGrowSynapseMutation <= growSynapseProbability) {

			const std::pair<id_t, id_t> neuronIds = genotype.findConnectableNeurons(m_numInputs);

			if (neuronIds.first < -m_numInputs && neuronIds.second < -m_numInputs) {
//				std::cout << "Failed to find connectable neurons!" << std::endl;
			} else {
				growSynapseMutations[neuronIds].push_back(genotypeId);
			}
		}

		// split synapse mutation
		if (rollForSplitSynapseMutation <= splitSynapseProbability) {
			const innov_t randomSynapseGeneId = genotype.findSplittableSynapse();

			const auto& it = genotype.getSynapseGenes().find(randomSynapseGeneId);
			if (it != genotype.getSynapseGenes().end()) {
//				std::cout << "	Splitting synapse gene " << randomSynapseGeneId << std::endl;
				synapseSplitMutations[randomSynapseGeneId].push_back(genotypeId);
			}
		}
	}

	// apply the split synapse mutations
	for (const auto& synapseSplitIt : synapseSplitMutations) {
		const innov_t synapseToSplitId = synapseSplitIt.first;
		const auto& genotypeIds = synapseSplitIt.second;

//		std::cout << "Splitting synapse " << synapseToSplitId << " in genotypes: ";

		for (const id_t genotypeId : genotypeIds) {
//			std::cout << genotypeId << " ";

			auto& genotype = m_genotypes[genotypeId];
			genotype.splitSynapse(m_innovationIndex+1, m_innovationIndex+2, m_neuronIndex, synapseToSplitId);
		}
		++m_neuronIndex;
		m_innovationIndex += 2;

//		std::cout << std::endl;
	}

	// apply the grow synapse mutations
	for (const auto& growSynapseIt : growSynapseMutations) {
		const std::pair<id_t, id_t>& neuronIds = growSynapseIt.first;
		const id_t startId = neuronIds.first;
		const id_t endId = neuronIds.second;
		const auto& genotypeIds = growSynapseIt.second;

//		std::cout << "Growing synapse " << startId << " -> " << endId << " in genotypes: ";

		++m_innovationIndex;
		for (const id_t genotypeId : genotypeIds) {
//			std::cout << genotypeId << " ";

			auto& genotype = m_genotypes[genotypeId];
			genotype.addSynapseGene(m_innovationIndex, Genotype::getRandomWeight(), startId, endId);
		}
//		std::cout << std::endl;
	}

	// mutate weights and biases
	for (auto& genotypeIt : m_genotypes) {
		const id_t genotypeId = genotypeIt.first;
		auto& genotype = genotypeIt.second;

		// champions should be unchanged
		if (m_championIds.count(genotypeId)) {
//			std::cout << "Champion " << genotypeId << " has its weights unchanged" << std::endl;
			continue;
		}

//		std::cout << "Weight mutation step for " << genotypeId << std::endl;
		double rollForMutation;

//		std::cout << "	Mutating synapses: ";
		for (const auto& synapseGeneIt : genotype.getSynapseGenes()) {
			rollForMutation =  NumberGenerator::getASym(1.0);

			if (rollForMutation <= mutateSynapseWeightProbability) {
				const innov_t synapseGeneId = synapseGeneIt.first;
				const auto& synapseGene = synapseGeneIt.second;

//				std::cout << synapseGeneId << " ";

				rollForMutation =  NumberGenerator::getASym(1.0);
				if (rollForMutation <= shiftSynapseWeightProbability) {
					const double oldWeight = synapseGene.getWeight();
					const double weight = oldWeight + NumberGenerator::getSym(0.1);
					genotype.setSynapseWeight(synapseGeneId, weight);
				} else {
					const double weight = Genotype::getRandomWeight();
					genotype.setSynapseWeight(synapseGeneId, weight);
				}
			}
		}
//		std::cout << std::endl;

//		std::cout << "Bias mutation step for " << genotypeId << std::endl;

//		std::cout << "	Mutating neurons: ";
		for (const auto &neuronGeneIt : genotype.getNeuronGenes()) {
			rollForMutation = NumberGenerator::getASym(1.0);

			if (rollForMutation <= mutateNeuronBiasProbability) {
				const id_t neuronGeneId = neuronGeneIt.first;
				const double oldBias = neuronGeneIt.second;

//				std::cout << neuronGeneId << " ";

				rollForMutation = NumberGenerator::getASym(1.0);
				if (rollForMutation <= shiftNeuronBiasProbability) {
					const double bias = oldBias + NumberGenerator::getSym(0.05);
					genotype.setNeuronBias(neuronGeneId, bias);
				} else {
					const double bias = Genotype::getRandomBias();
					genotype.setNeuronBias(neuronGeneId, bias);
				}
			}
		}
//		std::cout << std::endl;
	}

//	std::cout << "Mutation finished!" << std::endl << std::endl << std::endl;
}

void GenePool::constructNextGeneration() {
//	std::cout << "Constructing next generation!" << std::endl;

	speciate();

	removeWeakerGenotypes();

	removeWeakerSpecies();

	assignOffspringToSpecies();

	mateGenotypes();

	mutateGenotypes();

	for (auto& genotypeIt : m_genotypes) {
		auto& genotype = genotypeIt.second;
		genotype.orderNodes();
	}

	++m_generationId;
}
