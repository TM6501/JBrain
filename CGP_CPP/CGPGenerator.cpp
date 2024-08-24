#include "pch.h"
#include "CGPGenerator.h"
#include "CGPFunctions.h"
#include "CGPFFANNIndividual.h"
#include "StandardCGPIndividual.h"
#include "StandardCGPIndividual_impl.h"
#include "Enums.h"

#include "IAgent.h"

#include <vector>
#include <functional>
#include <map>
#include <iostream>
#include <tuple>

extern unsigned int DEBUG_LEVEL;

namespace CGP
{
	CGPGenerator::CGPGenerator() : m_initialized(false), m_agentType(CGP_TYPE::UNDEFINED)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;
	}

	// Inherited from Experiment::IAgentGenerator. Get YAML nodes corresponding
	// to the parameters that we'll need:
	bool CGPGenerator::loadConfigurationFromYaml(YAML::Node node)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Assume success until we find failure:
		bool configSuccess = true;		

		// Everything must be under "CGPIndividual":
		if (node["CGPIndividual"])
		{
			// Parameters are broken down to general and type-specific. Both
			// must be there:

			// Get the general parameters:
			if (node["CGPIndividual"]["CGPParameters"])
				m_cgpGeneralParameters = node["CGPIndividual"]["CGPParameters"];
			else
				configSuccess = false;

			// Only two types allowed for now:
			std::string cgpType = node["CGPType"].as<std::string>();
			if (cgpType == "CGPFFANN")
			{
				m_agentType = CGP_TYPE::CGPFFANN;
				m_versionSpecificParameters = node["CGPIndividual"]["CGPFFANNParameters"];
				m_nodeProcessingType = ""; // Not used with CGPFFANN
			}
			else if (cgpType == "StandardCGP")
			{
				m_agentType = CGP_TYPE::GENERAL;
				m_versionSpecificParameters = node["CGPIndividual"]["StandardCGP"];
				m_nodeProcessingType = node["CGPIndividual"]["ArgumentType"].as<std::string>();
			}
			else
			{
				std::cout << "CGP Type: " << cgpType << " not supported." << std::endl;
				configSuccess = false;
			}

			// Mutation parameters must be provided:
			if (node["CGPIndividual"]["MutationParameters"])
				m_mutationParameters = node["CGPIndividual"]["MutationParameters"];
			else
				configSuccess = false;
		}
		else
		{
			configSuccess = false;
		}

		if (configSuccess)
		{
			m_initialized = true;
		}

		return configSuccess;
	}

	Experiment::IAgent* CGPGenerator::getRandomCGPFFANNIndividual()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		Experiment::IAgent* agent = nullptr;

		std::vector<std::string> funcList;
		
		// Need to convert the YAML list to a C++ vector:
		for (YAML::const_iterator it = m_cgpGeneralParameters["functionList"].begin();
			it != m_cgpGeneralParameters["functionList"].end(); ++it)
		{
			funcList.push_back(it->as<std::string>());			
		}

		// Get our arguments ready for the actual constructor:
		bool constraintsUsed = false;
		double minConstraint = 0.0;
		double maxConstraint = 0.0;

		// If constraints are provided in the YAML, we use constraints:
		if (m_versionSpecificParameters["minConstraint"])
		{
			constraintsUsed = true;
			minConstraint = m_versionSpecificParameters["minConstraint"].as<double>();
			maxConstraint = m_versionSpecificParameters["maxConstraint"].as<double>();
		}

		// Create the individual:
		agent = new CGP::CGPFFANNIndividual(
			m_cgpGeneralParameters["inputSize"].as<int>(), // Input size
			m_cgpGeneralParameters["outputSize"].as<int>(), // Output size
			m_cgpGeneralParameters["rows"].as<int>(), // rows
			m_cgpGeneralParameters["columns"].as<int>(), // columns
			m_cgpGeneralParameters["columnsBack"].as<int>(),  // columns back for input        
			m_versionSpecificParameters["minWeight"].as<double>(),  // Min weight
			m_versionSpecificParameters["maxWeight"].as<double>(),  // maxWeight
			m_versionSpecificParameters["minNeuronInputCount"].as<int>(),  // Min neuron input count
			m_versionSpecificParameters["maxNeuronInputCount"].as<int>(),  // Max neuron input count
			constraintsUsed,  // Don't use constraints
			minConstraint,  // Min Constraint (not used)
			maxConstraint,  // max constraint (not used)
			m_versionSpecificParameters["minPreBias"].as<double>(),  // Minimum bias
			m_versionSpecificParameters["maxPreBias"].as<double>(),  // Maximum bias
			m_versionSpecificParameters["minPostBias"].as<double>(),  // Minimum bias
			m_versionSpecificParameters["maxPostBias"].as<double>(),  // Maximum bias
			funcList);  // Activation function list

		static_cast<AbstractCGPIndividual*>(agent)->randomize();

		return agent;
	}

	Experiment::IAgent* CGPGenerator::getRandomStandardCGPIndividual_double()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		Experiment::IAgent* agent = nullptr;
		std::vector<std::string> funcList;

		// Need to convert the YAML list to a C++ vector:
		for (YAML::const_iterator it = m_cgpGeneralParameters["functionList"].begin();
			it != m_cgpGeneralParameters["functionList"].end(); ++it)
		{
			funcList.push_back(it->as<std::string>());
		}
		
		agent = new CGP::StandardCGPIndividual<double>(
			m_cgpGeneralParameters["inputSize"].as<int>(), // Input size
			m_cgpGeneralParameters["outputSize"].as<int>(), // Output size
			m_cgpGeneralParameters["rows"].as<int>(), // rows
			m_cgpGeneralParameters["columns"].as<int>(), // columns
			m_cgpGeneralParameters["columnsBack"].as<int>(),  // columns back for input
			m_versionSpecificParameters["pRelevant"].as<bool>(), // pRelevant
			m_versionSpecificParameters["minP"].as<double>(), // minP
			m_versionSpecificParameters["maxP"].as<double>(), // maxP
			m_versionSpecificParameters["useConstraint"].as<bool>(), // useConstraint
			m_versionSpecificParameters["minConstraint"].as<double>(), // minConstraint
			m_versionSpecificParameters["maxConstraint"].as<double>(), // maxConstraint
			m_versionSpecificParameters["useScale"].as<bool>(), // useScale
			m_versionSpecificParameters["minScale"].as<double>(), // minScale
			m_versionSpecificParameters["maxScale"].as<double>(), // maxScale
			FUNCTION_TYPE::DOUBLE3_DOUBLE,
			funcList
		);

		static_cast<AbstractCGPIndividual*>(agent)->randomize();
		return agent;
	}

	Experiment::IAgent* CGPGenerator::getRandomStandardCGPIndividual_bool()
	{
		Experiment::IAgent* agent = nullptr;
		std::vector<std::string> funcList;

		// Need to convert the YAML list to a C++ vector:
		for (YAML::const_iterator it = m_cgpGeneralParameters["functionList"].begin();
			it != m_cgpGeneralParameters["functionList"].end(); ++it)
		{
			funcList.push_back(it->as<std::string>());
		}

		agent = new CGP::StandardCGPIndividual<bool>(
			m_cgpGeneralParameters["inputSize"].as<int>(), // Input size
			m_cgpGeneralParameters["outputSize"].as<int>(), // Output size
			m_cgpGeneralParameters["rows"].as<int>(), // rows
			m_cgpGeneralParameters["columns"].as<int>(), // columns
			m_cgpGeneralParameters["columnsBack"].as<int>(),  // columns back for input
			m_versionSpecificParameters["pRelevant"].as<bool>(), // pRelevant
			m_versionSpecificParameters["minP"].as<double>(), // minP
			m_versionSpecificParameters["maxP"].as<double>(), // maxP
			m_versionSpecificParameters["useConstraint"].as<bool>(), // useConstraint
			m_versionSpecificParameters["minConstraint"].as<double>(), // minConstraint
			m_versionSpecificParameters["maxConstraint"].as<double>(), // maxConstraint
			m_versionSpecificParameters["useScale"].as<bool>(), // useScale
			m_versionSpecificParameters["minScale"].as<double>(), // minScale
			m_versionSpecificParameters["maxScale"].as<double>(), // maxScale
			FUNCTION_TYPE::BOOL3_BOOL,
			funcList
			);

		static_cast<AbstractCGPIndividual*>(agent)->randomize();
		return agent;
	}

	Experiment::IAgent* CGPGenerator::getRandomIndividual()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		Experiment::IAgent* agent = nullptr;

		if (m_agentType == CGP_TYPE::GENERAL)
		{
			// Two types of arguments allowed for now:
			if (m_nodeProcessingType == "Double")
			{
				return getRandomStandardCGPIndividual_double();
			}
			else if (m_nodeProcessingType == "Bool")
			{
				return getRandomStandardCGPIndividual_bool();
			}
			else
			{
				std::cout << "Unrecognized node processing type: " << m_nodeProcessingType << std::endl;
				return nullptr;
			}			
		}

		else if (m_agentType == CGP_TYPE::CGPFFANN)
		{
			return getRandomCGPFFANNIndividual();			
		}

		// Don't need to worry about an else case. It gets stopped on initial
		// YAML-config read-in.
				
		return agent;		
	}

	std::vector<Experiment::IAgent*> CGPGenerator::getMutatedChildren(
		const std::vector<Experiment::IAgent*>& parents)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Get the mutation strategy:
		MUTATION_STRATEGY mutStrat = MUTATION_STRATEGY::UNDEFINED;

		if (m_mutationParameters["mutationType"].as<std::string>() == "PROBABILISTIC_PERGENE")
			mutStrat = MUTATION_STRATEGY::PROBABILISTIC_PERGENE;
		else if (m_mutationParameters["mutationType"].as<std::string>() == "PROBABILISTIC_PERVALUE")
			mutStrat = MUTATION_STRATEGY::PROBABILISTIC_PERVALUE;
		else if (m_mutationParameters["mutationType"].as<std::string>() == "ACTIVE_GENE")
			mutStrat = MUTATION_STRATEGY::ACTIVE_GENE;
		else
		{
			std::cout << "Unrecognized mutation strategy: " << m_mutationParameters["mutationType"].as<std::string>() << std::endl;
			exit(-1);
		}

		std::unordered_map<std::string, double> mutParams;

		for (YAML::const_iterator it = m_mutationParameters["mutationParameters"].begin();
			it != m_mutationParameters["mutationParameters"].end(); ++it)
		{
			mutParams.insert({ it->first.as<std::string>(), it->second.as<double>() });
		}

		// CGP only does 1-parent to 1-child mutation:
		AbstractCGPIndividual* parent = static_cast<AbstractCGPIndividual*>(parents[0]);
		Experiment::IAgent* child = static_cast<Experiment::IAgent*>(
			parent->getOneMutatedChild(mutStrat, mutParams));
		
		std::vector<Experiment::IAgent*> retChildren{ child };
		return retChildren;
	}
}