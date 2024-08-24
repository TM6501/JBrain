#include "pch.h"
#include "JBrainCGPIndividual.h"
#include "BasicGene.h"
#include "CGPFunctions.h"
#include <random>
#include <iostream>
#include "json.hpp"
using json = nlohmann::json;

namespace CGP
{
	JBrainCGPIndividual::JBrainCGPIndividual(
		unsigned int inputSize,
		unsigned int outputSize,
		unsigned int rows,
		unsigned int columns,
		unsigned int colsBack,
		float minP, float maxP,		
		std::vector<std::string> functionStringList,
		std::vector<std::function<double(double, double, double)> > functionList,
		float minConstraint = -2.0, float maxConstraint = 2.0,
		bool useConstraint = true)
		: AbstractCGPIndividual(static_cast<int>(inputSize),
			static_cast<int>(outputSize),
			static_cast<int>(rows),
			static_cast<int>(columns),
			static_cast<int>(colsBack),
			-1), // We don't use colsForward
		m_functionStringList(functionStringList),
		m_functionList(functionList),
		m_genome(),
		m_activeGenes(),
		m_minP(minP),
		m_maxP(maxP),
		m_minConstraint(minConstraint),  // Default to using a [-2.0, 2.0] constraint:
		m_maxConstraint(maxConstraint),
		m_useConstraint(useConstraint)
	{
	}

	JBrainCGPIndividual::JBrainCGPIndividual(const JBrainCGPIndividual& other)
		: JBrainCGPIndividual(other.m_inputSize, other.m_outputSize, other.m_rows,
			other.m_columns, other.m_columnsBack, other.m_minP, other.m_maxP,
			other.m_functionStringList, other.m_functionList, other.m_minConstraint,
			other.m_maxConstraint, other.m_useConstraint)
	{
		m_genome = other.m_genome;
		m_activeGenes = other.m_activeGenes;
	}

	JBrainCGPIndividual* JBrainCGPIndividual::getCGPIndividualFromJson(json& j)
	{
		// Confirm we're reading something good:
		if (!j["individualType"].is_string() ||
			j["individualType"].get<std::string>() != "JBrainCGPIndividual")
			return nullptr;

		// Now, assume everything is correct:
		unsigned int inputSize = j["inputSize"].get<unsigned int>();
		unsigned int outputSize = j["outputSize"].get<unsigned int>();
		unsigned int genomeSize = j["genoTypeSize"].get<unsigned int>();
		float minP = j["minP"].get<float>();
		float maxP = j["maxP"].get<float>();
		bool useConstraint = j["useConstraint"].get<bool>();
		float minConstraint = j["minConstraint"].get<float>();
		float maxConstraint = j["maxConstraint"].get<float>();

		// Build our function lists:
		std::vector<std::string> functionStringList = j["functionList"];
		std::vector<std::function<double(double, double, double)> > functionList;

		for (unsigned int i = 0; i < functionStringList.size(); ++i)
		{
			functionList.push_back(CGPFunctions::doubleIn_doubleOut::
				getFuncFromString(functionStringList[i]));
		}

		JBrainCGPIndividual* retVal = new JBrainCGPIndividual(
			inputSize,
			outputSize,
			1, // rows
			genomeSize, // columns
			genomeSize, // colsBack,
			minP,
			maxP,
			functionStringList,
			functionList,
			minConstraint,
			maxConstraint,
			useConstraint);

		// Make sure it is empty and reserve enough space:
		retVal->m_genome.clear();
		retVal->m_genome.reserve(genomeSize);

		// Read in each gene separately:
		GENE_TYPE gType;
		int x;
		int y;
		int func;
		double p;
		for (auto& elem : j["genes"])
		{
			x = -1;
			y = -1;
			func = -1;
			p = -1;
			gType = CGP::StringToGeneType(elem["type"].get<std::string>());
			
			// Input won't have anything else defined:
			if (gType != GENE_TYPE::INPUT)
			{ 
				x = elem["X"].get<int>();
			}

			// Processing has the rest:
			if (gType == GENE_TYPE::PROCESSING)
			{
				y = elem["Y"].get<int>();
				func = elem["F"].get<int>();
				p = elem["P"].get<double>();
			}

			retVal->m_genome.push_back(NoScaleGene(gType, x, y, func, p));
		}

		retVal->calculateActiveGenes();
		return retVal;
	}

	JBrainCGPIndividual::~JBrainCGPIndividual()
	{}

	void JBrainCGPIndividual::calculateActiveGenes()
	{
		// Start with a blank slate:
		m_activeGenes.clear();

		// A gene is active if it is used directly or indirectly by the
		// outputs from the CGP. Start with all of the outputs:
		for (int i = 0; i < m_outputSize; ++i)
		{
			m_activeGenes.push_back(int(m_genome.size() - i - 1));
		}

		// Go through all genes already added and add those they depend upon.
		// m_activeGenes will grow during this loop:
		std::vector<int> geneNumbsToAdd = {};

		for (unsigned int i = 0; i < m_activeGenes.size(); ++i)
		{
			auto gene = m_genome[m_activeGenes[i]];
			geneNumbsToAdd.clear();

			// Three gene types: Input, processing, and output:
			// Nothing to do with input genes:
			if (gene.type == GENE_TYPE::INPUT)
				continue;

			// X is used by output genes:
			else if (gene.type == GENE_TYPE::OUTPUT)
				geneNumbsToAdd.push_back(gene.X);

			// X and Y are used by processing genes:
			else if (gene.type == GENE_TYPE::PROCESSING)
			{
				geneNumbsToAdd.push_back(gene.X);
				geneNumbsToAdd.push_back(gene.Y);
			}

			else // Undefined gene type? Should never get here.
				assert(false);

			// For every gene used by this gene, if it hasn't already
			// been added to our list of active genes, do it now:
			for (unsigned int j = 0; j < geneNumbsToAdd.size(); ++j)
			{
				if (std::find(m_activeGenes.begin(), m_activeGenes.end(),
					geneNumbsToAdd[j]) == m_activeGenes.end())
				{
					m_activeGenes.push_back(geneNumbsToAdd[j]);
				}
			}
		}

		// Sort the active genes to make them easier to use elsewhere:
		std::sort(m_activeGenes.begin(), m_activeGenes.end());
	}

	void JBrainCGPIndividual::randomize()
	{
		// Cleare the genome:
		m_genome.clear();

		// Add all inputs:
		for (int i = 0; i < m_inputSize; ++i)
		{
			// Most parameters of input genes are irrelevant:
			m_genome.push_back(NoScaleGene(GENE_TYPE::INPUT, -1, -1, -1, 1.0));
		}

		// Create our random number generators:        
		static std::random_device rd;
		static std::mt19937 rng(rd());
		std::uniform_int_distribution<int> funcDistribution(0, int(m_functionList.size()) - 1);
		std::uniform_real_distribution<double> pDistribution(m_minP, m_maxP);

		for (int i = 0; i < (m_rows * m_columns); ++i)
		{
			auto inX = getValidInputNodeNumber(i + m_inputSize);
			auto inY = getValidInputNodeNumber(i + m_inputSize);
			auto inFuncNum = funcDistribution(rng);
			auto inP = pDistribution(rng);
			
			m_genome.push_back(NoScaleGene(GENE_TYPE::PROCESSING, inX,
				inY, inFuncNum, inP));
		}

		// Add the output genes:
		for (int i = 0; i < m_outputSize; ++i)
		{
			// X is all that matters for an output gene:
			auto inX = getValidInputNodeNumber(i + m_inputSize + (m_rows * m_columns));
			m_genome.push_back(NoScaleGene(GENE_TYPE::OUTPUT, inX, -1, -1, 1.0));
		}

		// Make sure we did everything right:
		assert(m_genome.size() == static_cast<unsigned int>(m_inputSize + m_outputSize + (m_rows * m_columns)));

		// Calculate the active genes. This needs to be done any time something
		// happened that could have changed them:
		calculateActiveGenes();
	}

	std::vector<double> JBrainCGPIndividual::calculateOutputs(const std::vector<double>& inputs)
	{
		// Track the outputs of each CGP node to use as inputs to others:
		std::unordered_map<int, double> nodeOutputs = {};

		// The outputs of input-type nodes are simply the values passed to us:
		for (unsigned int i = 0; i < inputs.size(); ++i)
			nodeOutputs.emplace(i, inputs[i]);

		// Step through our active genes, calculating each of their outputs.
		for (unsigned int i = 0; i < m_activeGenes.size(); ++i)
		{
			auto gene = m_genome[m_activeGenes[i]];

			// Ignore input genes:
			if (gene.type == GENE_TYPE::INPUT)
				continue;

			// Output genes just get the output from their X gene number:
			else if (gene.type == GENE_TYPE::OUTPUT)
			{
				// The value needs to be there:
				assert(nodeOutputs(gene.X) != nodeOutputs.end());
				nodeOutputs.emplace(m_activeGenes[i], nodeOutputs[gene.X]);
			}

			// Must be a processing node:
			else
			{
				// Input values should be there:
				assert(nodeOutputs.find(gene.X) != nodeOutputs.end());
				assert(nodeOutputs.find(gene.Y) != nodeOutputs.end());

				double inX = nodeOutputs[gene.X];
				double inY = nodeOutputs[gene.Y];
				double calculatedOutput;

				// Calculate and constrain the outputs:
				calculatedOutput = m_functionList[gene.F](inX, inY, gene.P);
				calculatedOutput = constrain(calculatedOutput);

				// Put the value into the node outputs:
				nodeOutputs.emplace(m_activeGenes[i], calculatedOutput);
			}
		}

		// Push all of the output values to the return vector:
		std::vector<double> retValue;

		for (unsigned int i = static_cast<unsigned int>(m_genome.size()) - m_outputSize; i < m_genome.size(); ++i)
		{
			retValue.push_back(nodeOutputs[i]);
		}

		return retValue;
	}

	double JBrainCGPIndividual::constrain(double value)
	{
		if (!m_useConstraint)
			return value;
		else  // Constrain the value between min and max:
			return std::min(m_maxConstraint, std::max(m_minConstraint, static_cast<float>(value)));
	}

	AbstractCGPIndividual* JBrainCGPIndividual::getOneMutatedChild(
		MUTATION_STRATEGY strategy,
		std::unordered_map<std::string, double> parameters)
	{
		// Create a new individual with our current parameters:
		AbstractCGPIndividual* retVal = new JBrainCGPIndividual(
			m_inputSize, m_outputSize, m_rows, m_columns,
			m_columnsBack, m_minP, m_maxP, m_functionStringList,
			m_functionList, m_minConstraint, m_maxConstraint,
			m_useConstraint);

		// Copy over parts of ourselves:
		dynamic_cast<JBrainCGPIndividual*>(retVal)->m_genome = m_genome;
		dynamic_cast<JBrainCGPIndividual*>(retVal)->m_activeGenes = m_activeGenes;
		
		// Tell the copy to mutate itself. Parameters is unused by this
		// variation on a CGP individual:
		dynamic_cast<JBrainCGPIndividual*>(retVal)->mutateSelf(strategy);

		// And return. The caller is responsible for deleting this pointer:
		return retVal;
	}

	void JBrainCGPIndividual::mutateSelf(MUTATION_STRATEGY strategy)
	{
		// Single gene active gene mutation is assumed and used here.
		static int activeGenesToModify = 1;
		int activeGenesModified = 0;

		// Create our random number generators:
		static std::random_device rd;
		static std::mt19937 rng(rd());
		std::uniform_int_distribution<int> funcDistribution(0, int(m_functionList.size()) - 1);
		std::uniform_real_distribution<double> pDistribution(m_minP, m_maxP);
		std::uniform_int_distribution<int> geneSelectionDistribution(
			m_inputSize, m_inputSize + m_outputSize + (m_rows * m_columns) - 1);
		
		// Build the list of gene parts we'll mutate:        
		std::vector<std::string> genePartList = { "X", "Y", "P" };

		// Unlikely, but possible that we only have one function for every gene:
		if (m_functionList.size() > 1)
			genePartList.push_back("F");

		// To choose randomly the part of the gene to modify:
		std::uniform_int_distribution<int> genePartDistribution(0,
			static_cast<int>(genePartList.size()) - 1);

		while (activeGenesModified < activeGenesToModify)
		{
			// Select the gene to modify and what part of that gene:
			int toModify = geneSelectionDistribution(rng);
			int genePart = genePartDistribution(rng);

			// Override: If we're modifying an output gene, we must modify X:
			if (m_genome[toModify].type == GENE_TYPE::OUTPUT)
				genePart = 0;

			// Mutate X:
			if (genePartList[genePart] == "X")
			{
				// To ensure the gene actually changed, keep
				// selecting new values until it changes. Exclude the extremely
				// rare case of having only a single input and being the first
				// processing node since this could create an infinite loop:
				if (!(m_inputSize < 2 && toModify == 1))
				{
					int startX = m_genome[toModify].X;
					while (startX == m_genome[toModify].X)
						m_genome[toModify].X = getValidInputNodeNumber(toModify);
				}
			}
			// Mutate Y, almost identical to X
			else if (genePartList[genePart] == "Y")
			{
				if (!(m_inputSize < 2 && toModify == 1))
				{
					int startY = m_genome[toModify].Y;
					while (startY == m_genome[toModify].Y)
						m_genome[toModify].Y = getValidInputNodeNumber(toModify);
				}
			}
			// Mutate the function:
			else if (genePartList[genePart] == "F")
			{
				int startFunc = m_genome[toModify].F;

				// Repeat until it changes:
				while (startFunc == m_genome[toModify].F)
					m_genome[toModify].F = funcDistribution(rng);
			}
			// Only other option: P
			else
			{
				m_genome[toModify].P = pDistribution(rng);
			}

			// If this was an active gene, increment our active gene count:
			if (std::find(m_activeGenes.begin(), m_activeGenes.end(), toModify)
				!= m_activeGenes.end())
			{
				++activeGenesModified;
			}
		}
	}

	void JBrainCGPIndividual::printGenotype()
	{
		int outputWidth = 2;
		if (m_genome.size() > 99)
			outputWidth = 3;
		if (m_genome.size() > 999)
			outputWidth = 4;

		if (m_genome.size() == 0)
		{
			std::cout << "No genes defined." << std::endl;
		}
		else
		{
			std::cout << std::setprecision(2) << getPercentageNodesUsed() 
				      << "% genes are active." << std::endl;

			for (unsigned int i = 0; i < m_genome.size(); ++i)
			{
				auto gene = m_genome[i];
				std::cout << std::setw(outputWidth) << i;

				// Input, no information needed other than type:
				if (gene.type == GENE_TYPE::INPUT)
					std::cout << " - Input" << std::endl;

				// Output, just need to know where it is getting its information:
				else if (gene.type == GENE_TYPE::OUTPUT)
					std::cout << " - Output from gene " << gene.X << std::endl;

				else // Processing. Need all information
				{
					std::cout << " -  X:" << std::setw(outputWidth) << gene.X;
					std::cout << "   Y:" << std::setw(outputWidth) << gene.Y;
					std::cout << "   P:" << std::setprecision(3) << gene.P;
					std::cout << "   F:" << m_functionStringList[gene.F];
					std::cout << std::endl;
				}
			}
		}

		std::cout << "Active genes: ";
		for (int i : m_activeGenes)
			std::cout << i << " ";
		std::cout << std::endl;
	}

	void JBrainCGPIndividual::resetForNewTimeSeries()
	{
		// Nothing needs to be done... ?
	}

	double JBrainCGPIndividual::getPercentageNodesUsed()
	{
		// Active genes are the only ones used:
		return (double(m_activeGenes.size()) / double(m_genome.size())) * 100.0;
	}

	void JBrainCGPIndividual::writeSelfToJson(json& j)
	{
		// Standard to all CGP:
		j["individualType"] = "JBrainCGPIndividual";
		j["inputSize"] = m_inputSize;
		j["outputSize"] = m_outputSize;
		j["genoTypeSize"] = m_genome.size();

		// P is always used, just provide min and max:
		j["minP"] = m_minP;
		j["maxP"] = m_maxP;

		// Constraint information:
		j["useConstraint"] = m_useConstraint;
		j["minConstraint"] = m_minConstraint;
		j["maxConstraint"] = m_maxConstraint;

		// Function list as strings. Function pointers don't carry over
		// from one execution to the next:
		j["functionList"] = json::array();
		for (unsigned int i = 0; i < m_functionStringList.size(); ++i)
			j["functionList"][i] = m_functionStringList[i];

		// Our genome, one gene at a time:
		j["genes"] = json::array();
		for (unsigned int i = 0; i < m_genome.size(); ++i)
		{
			// Gene type:
			j["genes"][i]["type"] = GeneTypeToString(m_genome[i].type);

			// All non-input need X:
			if (m_genome[i].type != GENE_TYPE::INPUT)
				j["genes"][i]["X"] = m_genome[i].X;

			// Processing needs everything (output handled by above check):
			if (m_genome[i].type == GENE_TYPE::PROCESSING)
			{
				j["genes"][i]["Y"] = m_genome[i].Y;
				j["genes"][i]["F"] = m_genome[i].F;
				j["genes"][i]["P"] = m_genome[i].P;
			}
		}
	}

	void JBrainCGPIndividual::readSelfFromJson(json& j)
	{
		// Not used.  getCGPIndividualFromJson is used instead.
	}

	void JBrainCGPIndividual::performOncePerEpochUpdates(
		std::vector<AbstractCGPIndividual*> population, double epochScore)
	{
		// Nothing to do.
	}

	JBrainCGPIndividual& JBrainCGPIndividual::operator=(const JBrainCGPIndividual& rhs)
	{
		// Prevent self-assignment:
		if (this == &rhs)
			return *this;

		// Copy over all variables:
		m_functionStringList = rhs.m_functionStringList;
		m_functionList = rhs.m_functionList;
		m_genome = rhs.m_genome;
		m_activeGenes = rhs.m_activeGenes;
		m_minP = rhs.m_minP;
		m_maxP = rhs.m_maxP;
		m_minConstraint = rhs.m_minConstraint;
		m_maxConstraint = rhs.m_maxConstraint;
		m_useConstraint = rhs.m_useConstraint;
	}

	bool JBrainCGPIndividual::operator==(const JBrainCGPIndividual& rhs)
	{
		if (this == &rhs)
			return true;
		
		// Apparently, if the == operator is defined for the elements,
		// C++ provides a vector<> == operator. We can compare them
		// as any other variable:
		return (m_functionStringList == rhs.m_functionStringList) &&
			(m_genome == rhs.m_genome) &&
			(m_activeGenes == rhs.m_activeGenes) &&
			(fabs(m_minP - rhs.m_minP) < FLT_EPSILON) &&
			(fabs(m_maxP - rhs.m_maxP) < FLT_EPSILON) &&
			(fabs(m_minConstraint - rhs.m_minConstraint) < FLT_EPSILON) &&
			(fabs(m_maxConstraint - rhs.m_maxConstraint) < FLT_EPSILON) &&
			(m_useConstraint == rhs.m_useConstraint);
	}

	// Comparison operators:
	// Component comparison operators:
	bool operator==(const NoScaleGene& lhs, const NoScaleGene& rhs)
	{
		if (&lhs == &rhs)
			return true;

		// Some parts are only required for specific gene types.
		// Compare only the required parts:
		if (lhs.type != rhs.type)
			return false;

		if (lhs.type == GENE_TYPE::INPUT || lhs.type == GENE_TYPE::UNDEFINED)
			return true;

		if (lhs.type == GENE_TYPE::OUTPUT)
			return lhs.X == rhs.X;

		// Processing gene, everything else required:
		return (lhs.Y == rhs.Y) &&
			(lhs.F == rhs.F) &&
			(fabs(lhs.P - rhs.P) < DBL_EPSILON);
	}

}  // End CGP namespace