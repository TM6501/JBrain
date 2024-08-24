#include "pch.h"
#include "CGPFFANNIndividual.h"
#include <random>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <string>

extern unsigned int DEBUG_LEVEL;

namespace CGP
{
	double CGPFFANNIndividual::BIAS_CHANGE_MINIMUM = 0.01;

	CGPFFANNIndividual::CGPFFANNIndividual(int inputSize, int outputSize,
		int rows, int columns, int colsBack, double minWeight,
		double maxWeight, int minNeuronInputCount, int maxNeuronInputCount,
		bool useConstraint, double minConstraint, double maxConstraint,
		double minPreBias, double maxPreBias, double minPostBias,
		double maxPostBias,
		std::vector<std::string> functionList)
		: AbstractCGPIndividual(inputSize, outputSize, rows, columns, colsBack, 0),
		m_minWeight(minWeight),
		m_maxWeight(maxWeight),
		m_minPreBias(minPreBias),
		m_maxPreBias(maxPreBias),
		m_mutatePreBias(m_maxPreBias - m_minPreBias >= BIAS_CHANGE_MINIMUM),
		m_minPostBias(minPostBias),
		m_maxPostBias(maxPostBias),
		m_mutatePostBias(m_maxPostBias - m_minPostBias >= BIAS_CHANGE_MINIMUM),
		m_useConstraint(useConstraint),
		m_minConstraint(minConstraint),
		m_maxConstraint(maxConstraint),
		m_minInputCount(minNeuronInputCount),
		m_maxInputCount(maxNeuronInputCount),
		m_functionStringList(functionList)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Fill in our actual functions with those represented by the input strings:
		for (unsigned int i = 0; i < m_functionStringList.size(); ++i)
		{
			m_functionList.push_back(
				CGPFunctions::neuronActivation::getFuncFromString(
					m_functionStringList[i]));
			m_activeFunctions.push_back(i);
		}
	}

	void CGPFFANNIndividual::randomize()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Clear all neurons:
		m_genotype.clear();

		// Add the input neurons:
		for (int i = 0; i < m_inputSize; ++i)
		{
			// Input neurons don't need any values defined:
			m_genotype.push_back(CGPNeuron(GENE_TYPE::INPUT));
		}

		// Create our random number generators:        
		static std::random_device rd;
		static std::mt19937 rng(rd());		
		std::uniform_int_distribution<int> inputCountDistribution(m_minInputCount, m_maxInputCount);
		std::uniform_real_distribution<double> weightDistribution(m_minWeight, m_maxWeight);
		std::uniform_real_distribution<double> preBiasDistribution(m_minPreBias, m_maxPreBias);
		std::uniform_real_distribution<double> postBiasDistribution(m_minPostBias, m_maxPostBias);
		std::uniform_int_distribution<int> funcDistribution(0, int(m_activeFunctions.size()) - 1);

		// Add all of the processing nodes:
		for (int i = 0; i < (m_rows * m_columns); ++i)
		{
			auto tempNeuron = CGPNeuron(GENE_TYPE::PROCESSING);

			int numberInputs = inputCountDistribution(rng);
			for (int j = 0; j < numberInputs; ++j)
			{
				// Every inputs requires an input node number and weight:
				tempNeuron.m_inputs.push_back(getValidInputNodeNumber(i + m_inputSize));
				tempNeuron.m_weights.push_back(weightDistribution(rng));
			}

			// Add the bias:
			tempNeuron.m_preBias = preBiasDistribution(rng);
			tempNeuron.m_postBias = postBiasDistribution(rng);

			// Add the function (often there will only be a single possibility):
			tempNeuron.m_function = m_activeFunctions[funcDistribution(rng)];

			m_genotype.push_back(tempNeuron);
		}

		// Add all output nodes:
		for (int i = 0; i < m_outputSize; ++i)
		{
			auto tempNeuron = CGPNeuron(GENE_TYPE::OUTPUT);

			// Output neurons just need a single input:
			tempNeuron.m_inputs.push_back(getValidInputNodeNumber(i + m_inputSize + (m_rows * m_columns)));
			
			m_genotype.push_back(tempNeuron);
		}

		// This should now be the size of our genome:
		assert(m_genotype.size() == static_cast<unsigned int> (m_inputSize + m_outputSize + (m_rows * m_columns)));

		// Always recalculate our active genes when we do something that
		// could change them:
		calculateActiveGenes();
	}

	AbstractCGPIndividual* CGPFFANNIndividual::getOneMutatedChild(
		MUTATION_STRATEGY strategy, std::unordered_map<std::string, double> parameters)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Copy ourselves:
		AbstractCGPIndividual* retVal = new CGPFFANNIndividual(
			m_inputSize, m_outputSize, m_rows, m_columns, m_columnsBack, m_minWeight,
			m_maxWeight, m_minInputCount, m_maxInputCount, m_useConstraint,
			m_minConstraint, m_maxConstraint, m_minPreBias, m_maxPreBias,
			m_minPostBias, m_maxPostBias, m_functionStringList);

		static_cast<CGPFFANNIndividual*>(retVal)->m_genotype = m_genotype;
		static_cast<CGPFFANNIndividual*>(retVal)->m_activeGenes = m_activeGenes;
		static_cast<CGPFFANNIndividual*>(retVal)->m_activeFunctions = m_activeFunctions;

		static_cast<CGPFFANNIndividual*>(retVal)->mutateSelf(strategy, parameters);

		return retVal;
	}

	void CGPFFANNIndividual::performOncePerEpochUpdates(
		std::vector<AbstractCGPIndividual*> population, double epochScore)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;
	}

	void CGPFFANNIndividual::writeSelfToJson(json& j)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		j["individualType"] = "CGPFFANNIndividual";
		j["inputSize"] = m_inputSize;
		j["outputSize"] = m_outputSize;
		j["genoTypeSize"] = m_genotype.size();
		j["functionList"] = json::array();
		
		// Add the list of functions:
		for (unsigned int i = 0; i < m_functionStringList.size(); ++i)
			j["functionList"][i] = m_functionStringList[i];

		j["activeFunctionList"] = json::array();
		for (unsigned int i = 0; i < m_activeFunctions.size(); ++i)
			j["activeFunctionList"][i] = m_activeFunctions[i];

		j["genes"] = json::array();
		for (unsigned int i = 0; i < m_genotype.size(); ++i)
		{
			j["genes"][i]["type"] = static_cast<int>(m_genotype[i].m_geneType);
			j["genes"][i]["preBias"] = m_genotype[i].m_preBias;
			j["genes"][i]["postBias"] = m_genotype[i].m_postBias;
			j["genes"][i]["function"] = m_genotype[i].m_function;
			j["genes"][i]["inputs"] = json::array();
			j["genes"][i]["weights"] = json::array();

			for (unsigned int k = 0; k < m_genotype[i].m_inputs.size(); ++k)
			{
				j["genes"][i]["inputs"][k] = m_genotype[i].m_inputs[k];
			}

			for (unsigned int k = 0; k < m_genotype[i].m_weights.size(); ++k)
			{
				j["genes"][i]["weights"][k] = m_genotype[i].m_weights[k];
			}
		}
	}

	void CGPFFANNIndividual::readSelfFromJson(json& j)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		std::string cgpType = j["individualType"];
		if (cgpType != "CGPFFANNIndividual")
		{
			std::cout << "Wrong type! Can't read in " << cgpType << std::endl;
			return;
		}
		m_inputSize = j["inputSize"];
		m_outputSize = j["outputSize"];
		
		// Get the list of functions:
		m_functionStringList.clear();
		for (unsigned int i = 0; i < j["functionList"].size(); ++i)
			m_functionStringList.push_back(j["functionList"][i].get<std::string>());

		// Get the list of active functions:
		m_activeFunctions.clear();
		for (unsigned int i = 0; i < j["activeFunctionList"].size(); ++i)
			m_activeFunctions.push_back(j["activeFunctionList"][i].get<unsigned int>());

		int genotypeSize = j["genoTypeSize"];
		m_genotype.clear();
		for (int i = 0; i < genotypeSize; ++i)
		{	
			GENE_TYPE neuronType = static_cast<GENE_TYPE>(j["genes"][i]["type"]);
			CGPNeuron tempNeuron(neuronType);
			if (neuronType != GENE_TYPE::INPUT)
			{
				tempNeuron.m_inputs.clear();
				for (unsigned int k = 0; k < j["genes"][i]["inputs"].size(); ++k)
					tempNeuron.m_inputs.push_back(j["genes"][i]["inputs"][k]);

				// tempNeuron.m_inputs = static_cast<std::vector<int>>(j["genes"][i]["inputs"]);
				if (neuronType != GENE_TYPE::OUTPUT)
				{
					tempNeuron.m_weights.clear();
					for (unsigned int k = 0; k < j["genes"][i]["weights"].size(); ++k)
						tempNeuron.m_weights.push_back(j["genes"][i]["weights"][k]);
					
					tempNeuron.m_function = j["genes"][i]["function"];
					tempNeuron.m_preBias = j["genes"][i]["preBias"];
					tempNeuron.m_postBias = j["genes"][i]["postBias"];
				}
			}
			
			m_genotype.push_back(tempNeuron);
		}

		calculateActiveGenes();
	}

	void CGPFFANNIndividual::printGenotype()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		for (int i = 0; i < static_cast<int>(m_genotype.size()); ++i)
		{
			auto gene = m_genotype[i];
			if (gene.m_geneType == GENE_TYPE::INPUT)
			{
				std::cout << i << ": INPUT" << std::endl;
			}

			else if (gene.m_geneType == GENE_TYPE::OUTPUT)
			{
				std::cout << i << ": OUTPUT from node " << gene.m_inputs[0] << std::endl;
			}

			else if (gene.m_geneType == GENE_TYPE::PROCESSING)
			{
				std::cout << i << ": PROCESSING " << int(gene.m_inputs.size())
					      << " inputs with the " << m_functionStringList[gene.m_function] 
					      << " function:" << std::endl;
				std::cout << "\tPreBias: " << std::setprecision(6) << gene.m_preBias << std::endl;
				std::cout << "\tPostBias: " << std::setprecision(6) << gene.m_postBias << std::endl;
				for (int j = 0; j < static_cast<int>(gene.m_inputs.size()); ++j)
				{
					std::cout << "\tInput " << std::setw(5) << gene.m_inputs[j] << " weight: "
						      << std::setprecision(6) << gene.m_weights[j] << std::endl;
				}
			}
			else
			{
				assert(false);  // Unrecognized Gene type
			}
		}
	}

	void CGPFFANNIndividual::resetForNewTimeSeries()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;
	}

	double CGPFFANNIndividual::getPercentageNodesUsed()
	{
		return (double(m_activeGenes.size()) / double(m_genotype.size())) * 100.0;
	}

	void CGPFFANNIndividual::calculateActiveGenes()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Sets don't have the [] operator, so we'll fake some of the 
		// set's functionality as we ensure values aren't double-added:
		m_activeGenes.clear();

		// All output genes are considered active:
		for (int i = 1; i <= m_outputSize; ++i)
		{
			m_activeGenes.push_back(static_cast<int>(m_genotype.size() - i));
		}

		// Go through all genes already added and add those they require.
		// m_activeGenes will grow during this loop:
		std::vector<int> geneNumsToAdd = {};
		for (unsigned int i = 0; i < m_activeGenes.size(); ++i)
		{
			auto gene = m_genotype[m_activeGenes[i]];
			geneNumsToAdd.clear();

			// Nothing to do with input genes:
			if (gene.m_geneType == GENE_TYPE::INPUT)
				continue;

			// Only need a single input from output genes:
			else if (gene.m_geneType == GENE_TYPE::OUTPUT)
			{
				geneNumsToAdd.push_back(gene.m_inputs[0]);
			}

			else // Processing, add all inputs:
			{
				for (unsigned int j = 0; j < gene.m_inputs.size(); ++j)
				{
					geneNumsToAdd.push_back(gene.m_inputs[j]);
				}
			}

			// For every gene to add, if it isn't in the active genes already,
			// add it:
			for (unsigned int j = 0; j < geneNumsToAdd.size(); ++j)
			{
				if (std::find(m_activeGenes.begin(), m_activeGenes.end(),
					geneNumsToAdd[j]) == m_activeGenes.end())
				{
					m_activeGenes.push_back(geneNumsToAdd[j]);
				}
			}
		}

		// Sort the list of active genes now that we have them all:
		std::sort(m_activeGenes.begin(), m_activeGenes.end());
	}
	
	void CGPFFANNIndividual::mutateSelf(MUTATION_STRATEGY strategy,
		std::unordered_map<std::string, double> parameters)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		if (strategy == MUTATION_STRATEGY::PROBABILISTIC_PERGENE)
			mutateSelf_probabilisticPerGene(parameters);
		else if (strategy == MUTATION_STRATEGY::PROBABILISTIC_PERVALUE)
			mutateSelf_probabilisticPerValue(parameters);
		else if (strategy == MUTATION_STRATEGY::ACTIVE_GENE)
			mutateSelf_activeGene(parameters);
		else
		{
			// Shouldn't reach here:
			assert(false);
		}

		calculateActiveGenes();
	}

	void CGPFFANNIndividual::mutateSelf_probabilisticPerValue(
		std::unordered_map<std::string, double> parameters)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Make sure the values we expect are there:
		assert(parameters.find("percentage") != parameters.end());

		// Get the percentage of values that we will modify:
		double percentage = parameters.find("percentage")->second;

		// Create our random number generators:
		static std::random_device rd;
		static std::mt19937 rng(rd());		
		std::uniform_real_distribution<double> weightDistribution(m_minWeight, m_maxWeight);
		std::uniform_real_distribution<double> preBiasDistribution(m_minPreBias, m_maxPreBias);
		std::uniform_real_distribution<double> postBiasDistribution(m_minPostBias, m_maxPostBias);
		std::uniform_real_distribution<double> mutateChanceDistribution(0.0, 100.0);

		// If there's only 1 function, we can't mutate it:
		std::uniform_int_distribution<int> funcDistribution(0, int(m_activeFunctions.size()) - 1);
		bool mutateFunction = (m_functionList.size() > 1);

		// Step through every value of every gene. Mutate it if the random chance dictates
		// that we should. Input genes can't be mutated:
		for (unsigned int i = m_inputSize; i < m_genotype.size(); ++i)
		{
			// Output genes only have a single input to mutate:
			if (m_genotype[i].m_geneType == GENE_TYPE::OUTPUT)
			{
				if (mutateChanceDistribution(rng) < percentage)
				{
					m_genotype[i].m_inputs[0] = getValidInputNodeNumber(i);
				}
			}

			// Processing genes. Check every value independently:
			else
			{
				// Each input and weight:
				for (int j = 0; j < static_cast<int>(m_genotype[i].m_inputs.size()); ++j)
				{
					// Input:
					if (mutateChanceDistribution(rng) < percentage)
					{
						// Keep selecting new inputs until it changes.
						// Exclude the extremely rare case of having only a single
						// input and being the first processing node. This
						// will create an infinite loop:					
						if (!(m_inputSize < 2 && i == 1))
						{
							int startInput = m_genotype[i].m_inputs[j];
							while (startInput == m_genotype[i].m_inputs[j])
								m_genotype[i].m_inputs[j] = getValidInputNodeNumber(i);
						}
					}

					// Weight:
					if (mutateChanceDistribution(rng) < percentage)
					{
						m_genotype[i].m_weights[j] = weightDistribution(rng);
					}
				}

				// Bias:
				if (m_mutatePreBias && mutateChanceDistribution(rng) < percentage)
				{
					m_genotype[i].m_preBias = preBiasDistribution(rng);
				}
				
				if (m_mutatePostBias && mutateChanceDistribution(rng) < percentage)
				{
					m_genotype[i].m_postBias = postBiasDistribution(rng);
				}

				// Function:
				if (mutateFunction && mutateChanceDistribution(rng) < percentage)
				{
					int startFunc = m_genotype[i].m_function;

					// Keep selecting until it changes:
					while (startFunc == m_genotype[i].m_function)
					{
						m_genotype[i].m_function = m_activeFunctions[funcDistribution(rng)];
					}
				}

				// Add input if we can:
				if (m_genotype[i].m_inputs.size() < m_maxInputCount)
				{
					if (mutateChanceDistribution(rng) < percentage)
					{
						m_genotype[i].m_inputs.push_back(getValidInputNodeNumber(i));
						m_genotype[i].m_weights.push_back(weightDistribution(rng));
					}
				}

				// Remove input in we can:
				if (m_genotype[i].m_inputs.size() > m_minInputCount)
				{
					if (mutateChanceDistribution(rng) < percentage)
					{
						// Decide on the input to remove:
						std::uniform_int_distribution<int> inputDistribution(0, int(m_genotype[i].m_inputs.size() - 1));
						int toRemove = inputDistribution(rng);

						// Erase it from the input node vector and weight vector:
						m_genotype[i].m_inputs.erase(m_genotype[i].m_inputs.begin() + toRemove);
						m_genotype[i].m_weights.erase(m_genotype[i].m_weights.begin() + toRemove);
					}
				}
			}
		}
	}

	void CGPFFANNIndividual::mutateSelf_probabilisticPerGene(
		std::unordered_map<std::string, double> parameters)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Still needs definition
	}

	void CGPFFANNIndividual::mutateSelf_activeGene(
		std::unordered_map<std::string, double> parameters)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Make sure the values we expect are there:
		assert(parameters.find("geneCount") != parameters.end());

		// Get the percentage of values that we will modify:
		int genesToModify = static_cast<int>(parameters.find("geneCount")->second + 0.49);
		int genesModified = 0;

		// Create our random number generators:
		// Create our random number generators:
		static std::random_device rd;
		static std::mt19937 rng(rd());		
		std::uniform_real_distribution<double> preBiasDistribution(m_minPreBias, m_maxPreBias);
		std::uniform_real_distribution<double> postBiasDistribution(m_minPostBias, m_maxPostBias);
		std::uniform_int_distribution<int> geneSelectionDistribution(
			m_inputSize, m_inputSize + m_outputSize + (m_rows * m_columns) - 1);
		std::uniform_real_distribution<double> weightDistribution(m_minWeight, m_maxWeight);

		// If there's only 1 function, we can't mutate it:
		std::uniform_int_distribution<int> funcDistribution(0, int(m_activeFunctions.size()) - 1);
		bool mutateFunction = (m_functionList.size() > 1);

		// Build the starting list of gene parts we can mutate:
		// Start with just pre-bias and post-bias:
		std::vector<std::string> startGenePartList = {};
		if (m_mutatePreBias)
			startGenePartList.push_back("PreB");
		
		if (m_mutatePostBias)
			startGenePartList.push_back("PostB");

		if (mutateFunction)
		{
			startGenePartList.push_back("F");  // Add function if it is available
		}
		
		while (genesModified < genesToModify)
		{
			// Select the gene to modify:
			int toModify = geneSelectionDistribution(rng);

			// If this was an active gene, increment how many active genes
			// we've modified:
			if (std::find(m_activeGenes.begin(), m_activeGenes.end(), toModify)
				!= m_activeGenes.end())
			{
				++genesModified;
			}

			// Output genes have only a single value to modify:
			if (m_genotype[toModify].m_geneType == GENE_TYPE::OUTPUT)
			{
				// Keep changing it until we get something different:
				auto startInput = m_genotype[toModify].m_inputs[0];
				while (startInput == m_genotype[toModify].m_inputs[0])
				{
					m_genotype[toModify].m_inputs[0] = getValidInputNodeNumber(toModify);
				}
			}

			else
			{
				// Build the list of parts of that gene that can be modified:
				std::vector<std::string> genePartList = startGenePartList;
				if (m_genotype[toModify].m_inputs.size() > m_minInputCount)
					genePartList.push_back("R");  // Add remove if it is an option.

				if (m_genotype[toModify].m_inputs.size() < m_maxInputCount)
					genePartList.push_back("A");  // Add add if it is an option.

				// Add every input and weight:
				for (unsigned int i = 0; i < m_genotype[toModify].m_inputs.size(); ++i)
				{
					genePartList.push_back("W" + std::to_string(i));
					genePartList.push_back("I" + std::to_string(i));
				}

				// Select randomly the gene part to modify:
				std::uniform_int_distribution<int> partDist(
					0, static_cast<int>(genePartList.size()) - 1);

				std::string mod = genePartList[partDist(rng)];

				if (mod == "PreB")
				{
					// Change the bias:
					m_genotype[toModify].m_preBias = preBiasDistribution(rng);
				}

				else if (mod == "PostB")
				{
					// Change the bias:
					m_genotype[toModify].m_postBias = postBiasDistribution(rng);
				}

				else if (mod == "F")
				{
					int startFunc = m_genotype[toModify].m_function;
					while (startFunc == m_genotype[toModify].m_function)
					{
						m_genotype[toModify].m_function = m_activeFunctions[funcDistribution(rng)];
					}
				}

				else if (mod == "R")
				{
					// Decide which input to remove:
					std::uniform_int_distribution<int> tempDist(0,
						static_cast<int>(m_genotype[toModify].m_inputs.size()) - 1);

					auto toRemove = tempDist(rng);

					// Remove the corresponding weight and input:
					m_genotype[toModify].m_inputs.erase(
						m_genotype[toModify].m_inputs.begin() + toRemove);

					m_genotype[toModify].m_weights.erase(
						m_genotype[toModify].m_weights.begin() + toRemove);
				}

				else if (mod == "A")
				{
					// Add an input:
					m_genotype[toModify].m_inputs.push_back(getValidInputNodeNumber(toModify));
					m_genotype[toModify].m_weights.push_back(weightDistribution(rng));
				}

				else if (mod[0] == 'W')
				{
					// Check which weight to modify:
					int toMod = std::stoi(mod.substr(1, mod.length()));
					m_genotype[toModify].m_weights[toMod] = weightDistribution(rng);
				}

				else if (mod[0] == 'I')
				{
					// Check which input to modify:
					int toMod = std::stoi(mod.substr(1, mod.length()));

					// Keep selecting new inputs until it changes.
					// Exclude the extremely rare case of having only a single
					// input and being the first processing node, which could create
					// an infinite loop:
					if (!(m_inputSize < 2 && toModify == 1))
					{
						int startIn = m_genotype[toModify].m_inputs[toMod];
						while (startIn == m_genotype[toModify].m_inputs[toMod])
							m_genotype[toModify].m_inputs[toMod] =
							getValidInputNodeNumber(toModify);
					}
				}
				else
				{
					// Shouldn't get here:
					assert(false);
				}
			}  // End else-must-be-a-processing-node
		}  // End while genesMod < genesToMod
	} // End function

	std::vector<double> CGPFFANNIndividual::calculateOutputs(const std::vector<double>& inputs)
	{
		// This function happens too much for standard debugging to be useful:
		if (DEBUG_LEVEL > 5)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// The outputs calculated by each node:
		std::unordered_map<int, double> nodeOutputs = {};

		for (unsigned int i = 0; i < inputs.size(); ++i)
		{
			nodeOutputs.emplace(i, inputs[i]);
		}

		// Calculate the output from each node as needed. Generally, most aren't needed:
		for (unsigned int i = 0; i < m_activeGenes.size(); ++i)
		{
			auto gene = m_genotype[m_activeGenes[i]];

			// Ignore input genes:
			if (gene.m_geneType == GENE_TYPE::INPUT)
			{
				continue;
			}

			// Grab the only input value:
			else if (gene.m_geneType == GENE_TYPE::OUTPUT)
			{
				// The value should be there:
				assert(nodeOutputs.find(gene.m_inputs[0]) != nodeOutputs.end());

				nodeOutputs.emplace(m_activeGenes[i], nodeOutputs[gene.m_inputs[0]]);
			}

			// Otherwise, do the processing:
			else
			{
				// Calculate the total value to pass into the activation function.
				// Pre-bias is part of that total that gets passed into the function:
				double weightedTotal = gene.m_preBias;

				for (unsigned int j = 0; j < gene.m_inputs.size(); ++j)
				{
					weightedTotal += nodeOutputs[gene.m_inputs[j]] * gene.m_weights[j];
				}

				// Get the output from the activation function:				
				double output = m_functionList[gene.m_function](weightedTotal);

				// Post-bias gets added after the calculation:
				weightedTotal += gene.m_postBias;

				// Constraint if needed:
				if (m_useConstraint)
				{
					output = constrain(output);
				}

				// Put the value in our outputs:
				nodeOutputs.emplace(m_activeGenes[i], output);
			}
		}

		// Push all of the output values to our return vector:
		std::vector<double> retValue;

		for (auto i = m_genotype.size() - m_outputSize; i < m_genotype.size(); ++i)
		{
			retValue.push_back(nodeOutputs[static_cast<int>(i)]);
		}

		return retValue;
	}

};  // End namespace CGP
