#include "Population.hpp"

#include <cmath>
#include <thread>

#include "Phenotype.hpp"
#include "util/NumberGenerator.hpp"


#include "Interface.hpp"
#include "mvc/Model.hpp"

#include <iostream>
#include <iomanip>
#include "rendering/DAGRenderer.hpp"

Population::Population(const uint64_t numGenotypes, const uint16_t numInputs, const uint16_t numOutputs)
	: m_genePool(numGenotypes, numInputs, numOutputs)
	, m_scoringVisualizer({0.4, 0.0, 0.6, 1.0})
	, m_speciesRecordViewer(0.0, 0.5, 0.4, 0.5)
	, m_fitnessRecordViewer(0.0, 0.0, 0.4, 0.5)
	, m_generationId(0)
	, m_topGenerationFitness(0)
	, m_fitnessRecord()
	, m_bestGenotypeIdInGeneration(-1) {
}

//double Population::scorePhenotype(const std::unique_ptr<Phenotype>& phenotypePtr) {
//	std::vector<std::pair<int, int>> inputVariants = {std::make_pair(0,0), std::make_pair(0,1), std::make_pair(1,0), std::make_pair(1,1)};
//	std::vector<double> desiredOutputs = {0.0, 1.0, 1.0, 0.0};
//
//	double totalError = 0;
////	std::cout << "Outputs: ";
//	for (size_t i(0); i < inputVariants.size(); ++i) {
//		const auto& inputVariant = inputVariants[i];
//		const double desiredOutput = desiredOutputs[i];
//
//		std::unordered_map<id_t, double> inputs;
//		inputs[-1] = inputVariant.first; inputs[-2] = inputVariant.second;// inputs[-3] = 1;
//		const auto& outputs = phenotypePtr->execute(inputs);
//
//		const double output = Phenotype::activateSigmoid(outputs.find(0)->second);
//		const double d = (desiredOutput - output);
//		const double error = d*d;
//		totalError += error;
//
////		std::cout << output << " -> " << error << " (" << inputVariant.first << ", " << inputVariant.second << ") ";
//	}
////	std::cout << std::endl;
////	std::cout << "	Total error: " << totalError << std::endl;
//
//	const double score = std::max(4.0 - totalError, 0.0);
//
//	return score;
//}

//double Population::scorePhenotype(const std::unique_ptr<Phenotype>& phenotypePtr) {
//    Model model(Interface::getRandomSnake());
//    uint64_t steps = 0;
//
//    while (!model.gameIsOver() && steps < 1000) {
//		const std::unordered_map<id_t, double> inputs = Interface::generateInputs(model);
//		const auto& outputs = phenotypePtr->execute(inputs);
//		Interface::interpretOutputs(model, outputs);
//
//		++steps;
//    }
//
//    return model.getScore();
//}

//double Population::scorePhenotype(const std::unique_ptr<Phenotype>& phenotypePtr) {
//    FBModel model;
//
//    while (!model.gameIsOver()) {
//		const std::unordered_map<id_t, double> inputs = FBInterface::generateInputs(model);
//		const auto& outputs = phenotypePtr->execute(inputs);
//		FBInterface::interpretOutputs(model, outputs);
//
//		if (model.getScore() > 1000) {
//			break;
//		}
//    }
//
//    return model.getScore();
//}

double Population::scorePhenotype(const std::unique_ptr<Phenotype>& phenotypePtr) {
//	std::cout << "	Scoring Phenotype!!" << std::endl;
    Model model;
    uint64_t steps = 0;

    while (!model.gameIsOver() && steps < 50000) {
		const std::unordered_map<id_t, double> inputs = Interface::generateInputs(model);
		const auto& outputs = phenotypePtr->execute(inputs);
		Interface::interpretOutputs(model, outputs);

		++steps;
//		std::cout << model.getScore() << std::endl;
    }

    return model.getScore();
}

void Population::scorePhenotypeShell(std::unordered_map<id_t, std::unique_ptr<Phenotype>>::iterator phenotypeIt, const std::unordered_map<id_t, size_t>& vectorIndices, std::vector<double>& genotypeScoresVector, const uint32_t iterationsThisThread) {
    for (uint32_t i(0); i < iterationsThisThread; ++i) {
    	const id_t genotypeId = phenotypeIt->first;
    	const std::unique_ptr<Phenotype>& phenotypePtr = phenotypeIt->second;

    	const size_t vectorIndex = vectorIndices.find(genotypeId)->second;
    	double score = scorePhenotype(phenotypePtr);
		genotypeScoresVector[vectorIndex] = score;
		++phenotypeIt;
    }
}

void Population::scoreGenerationThreaded() {
	// reset
	m_topGenerationFitness = 0;

	// find number of available threads
    const uint32_t numThreads = std::thread::hardware_concurrency();

    // create phenotypes
	auto& genotypes = m_genePool.getGenotypes();
	std::unordered_map<id_t, std::unique_ptr<Phenotype>> phenotypes;

	for (auto& genotypeIt : genotypes) {
		const id_t genotypeId = genotypeIt.first;
		auto& genotype = genotypeIt.second;

		std::unique_ptr<Phenotype> phenotypePtr = std::make_unique<Phenotype>(genotype);
		phenotypes.emplace(genotypeId, std::move(phenotypePtr));
	}

	// find number of phenotypes to be scored by each thread
    const uint32_t numPhenotypes = phenotypes.size();
    const uint32_t rawPhenotypesPerThread = numPhenotypes / numThreads;
    const uint32_t rawSum = numThreads * rawPhenotypesPerThread;
    uint32_t difference = numPhenotypes - rawSum;

    // create a conversion from indices to genotypeIds
	std::unordered_map<id_t, size_t> vectorIndices(numPhenotypes);
	std::unordered_map<id_t, Genotype>::iterator genotypeIt = genotypes.begin();
	for (size_t i(0); i < numPhenotypes; ++i) {
		const id_t genotypeId = genotypeIt->first;
		vectorIndices[genotypeId] = i;

		++genotypeIt;
	}

//	for (const auto& it : vectorIndices) {
//		std::cout << it.first << " -> " << it.second << std::endl;
//	}

	// setup threads
	std::unordered_map<id_t, std::unique_ptr<Phenotype>>::iterator phenotypeIt = phenotypes.begin();
	std::vector<double> genotypeScoresVector(numPhenotypes);
	std::vector<std::thread> threads;

	// start threads
    for (uint32_t i = 0; i < numThreads; ++i) {
    	uint32_t iterationsThisThread = rawPhenotypesPerThread;
    	if (difference > 0) {
    		++iterationsThisThread;
    		--difference;
    	}

        threads.emplace_back(Population::scorePhenotypeShell, phenotypeIt, std::cref(vectorIndices), std::ref(genotypeScoresVector), iterationsThisThread);
        for (uint32_t i(0); i < iterationsThisThread; ++i) {
        	++phenotypeIt;
        }
    }

    // join threads
    for (auto& thread : threads) {
        thread.join();
    }

    // extract scores from the vector and place them into an unordered_map
//	for (size_t i(0); i < genotypeScoresVector.size(); ++i) {
//		std::cout << i << " -> " << genotypeScoresVector[i] << std::endl;
//	}

    std::unordered_map<id_t, double> genotypeScores;
    for (const auto& it : vectorIndices) {
    	const id_t genotypeId = it.first;
    	const size_t vectorIndex = it.second;

    	genotypeScores[genotypeId] = genotypeScoresVector[vectorIndex];
    }

    // find best scores
	id_t bestGenotypeIdInGeneration = -1;
	for (const auto& scoreIt : genotypeScores) {
		const id_t genotypeId =	scoreIt.first;
		const double score = scoreIt.second;

		if (score > m_topGenerationFitness) {
			m_topGenerationFitness = score;
			bestGenotypeIdInGeneration = genotypeId;
		}

//		std::cout << "	Genotype " << genotypeId << " scored " << score << std::endl;

		genotypeScores[genotypeId] = score;
	}

	m_genePool.setGenotypeScores(genotypeScores);
	m_fitnessRecord.push_back(m_topGenerationFitness);

	int precision = std::numeric_limits<double>::max_digits10;
	std::cout << std::fixed << std::setprecision(precision) << "Top generation " << m_generationId << " fitness: " << m_topGenerationFitness << " / " << m_genePool.getGenotypeFitnessRecord() << " from genotype " << bestGenotypeIdInGeneration << std::endl;

	m_bestGenotypeIdInGeneration = bestGenotypeIdInGeneration;
}

void Population::scoreGeneration() {
	m_topGenerationFitness = 0;

	auto& genotypes = m_genePool.getGenotypes();

	std::unordered_map<id_t, std::unique_ptr<Phenotype>> phenotypes;
	for (auto& genotypeIt : genotypes) {
		const id_t genotypeId = genotypeIt.first;
		Genotype& genotype = genotypeIt.second;

		std::unique_ptr<Phenotype> phenotype = std::make_unique<Phenotype>(genotype);
		phenotypes.emplace(genotypeId, std::move(phenotype));
	}

	id_t bestGenotypeIdInGeneration = -1;
	std::unordered_map<id_t, double> genotypeScores;
	for (const auto& phenotypeIt : phenotypes) {
		const id_t genotypeId =	phenotypeIt.first;
		const std::unique_ptr<Phenotype>& phenotypePtr = phenotypeIt.second;

		const double score = scorePhenotype(phenotypePtr);

		if (score > m_topGenerationFitness) {
			m_topGenerationFitness = score;
			bestGenotypeIdInGeneration = genotypeId;
		}

//		std::cout << "	Genotype " << genotypeId << " scored " << score << std::endl;

		genotypeScores[genotypeId] = score;
	}

	m_genePool.setGenotypeScores(genotypeScores);
	m_fitnessRecord.push_back(m_topGenerationFitness);

	int precision = std::numeric_limits<double>::max_digits10;
	std::cout << std::fixed << std::setprecision(precision) << "Top generation " << m_generationId << " fitness: " << m_topGenerationFitness << " / " << m_genePool.getGenotypeFitnessRecord() << " from genotype " << bestGenotypeIdInGeneration << std::endl;

	m_bestGenotypeIdInGeneration = bestGenotypeIdInGeneration;
}

void Population::selection() {
	++m_generationId;
	m_genePool.constructNextGeneration();
}

void Population::perform(sf::RenderWindow& window, bool& showFitnessRecord, bool& showSpeciesRecord, bool& showGenotype, bool& showPerformance, double& msPerFrame) {
	if (m_genePool.getGenotypes().count(m_bestGenotypeIdInGeneration) == 0) {
		std::cerr << "There was no best genotype in this generation to perform!" << std::endl;
	}

	Genotype& bestGenotypeInGeneration = m_genePool.getGenotypes().find(m_bestGenotypeIdInGeneration)->second;
	m_scoringVisualizer.loadGenotype(window, bestGenotypeInGeneration);

	Phenotype phenotype(bestGenotypeInGeneration);

//	FBModel model;
//	Model model(Interface::getRandomSnake());
	Model model;
	uint64_t steps = 0;

	sf::Clock clock;

	bool skip = false;
	while (!model.gameIsOver() && window.isOpen()) {
		sf::Time elapsed = clock.getElapsedTime();
		if (elapsed.asMilliseconds() < msPerFrame) {
			continue;
		}
		clock.restart();

		sf::Event event;
		while (window.pollEvent(event)) {
			switch (event.type) {
			case sf::Event::Closed:
				window.close();
				break;
			case sf::Event::KeyPressed:
				if (event.key.code == sf::Keyboard::F) {
					showFitnessRecord = !showFitnessRecord;
				} else if (event.key.code == sf::Keyboard::S) {
					showSpeciesRecord = !showSpeciesRecord;
				} else if (event.key.code == sf::Keyboard::G) {
					showGenotype = !showGenotype;
				} else if (event.key.code == sf::Keyboard::P) {
					showPerformance = !showPerformance;
					skip = true;
				} else if (event.key.code == sf::Keyboard::Space) {
					skip = true;
				} else if (event.key.code == sf::Keyboard::A) {
					msPerFrame = std::max(1.0, msPerFrame - 1);
				} else if (event.key.code == sf::Keyboard::D) {
					msPerFrame = std::min(100.0, msPerFrame + 1);
				} else if (event.key.code == sf::Keyboard::Escape) {
					window.close();
				}
				break;
			default:
				break;
			}
		}

		if (skip) {
			window.clear();
			window.display();
			break;
		}

		const std::unordered_map<id_t, double> inputs = Interface::generateInputs(model);
//		const std::unordered_map<id_t, double> inputs = FBInterface::generateInputs(model);

		const auto& outputs = phenotype.execute(inputs);

		Interface::interpretOutputs(model, outputs);
//		FBInterface::interpretOutputs(model, outputs);

		++steps;

		window.clear();

		if (showGenotype) {
			m_scoringVisualizer.renderGenotype(window);
		}


		m_scoringVisualizer.renderScoreboard(window, model.getScore());
		m_scoringVisualizer.renderSimulation(window, model);

		if (showSpeciesRecord) {
			m_speciesRecordViewer.renderSpeciesRecord(window, m_genePool.getSpeciesRecord());
		}
		if (showFitnessRecord) {
			m_fitnessRecordViewer.renderFitnessRecord(window, m_fitnessRecord);
		}

		window.display();
	}
}
