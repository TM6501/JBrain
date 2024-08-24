#pragma once
#include "IAgentGenerator.h"
#include "Enums.h"
#include "yaml-cpp/yaml.h"

namespace CGP
{
	// The CGPGenerator class is here to satisify the requirements set out
	// by the full experiment for a class that can create individuals, mutate
	// them, and take all of its parameters from a YAML file.
	class CGPGenerator : public Experiment::IAgentGenerator
	{
	public:
		// Constructor does almost nothing. Parameters need to be read in via
		// the readInAllParameters function before the class is useful.
		CGPGenerator();

		// Inherited from Experiment::IAgentGenerator
		bool loadConfigurationFromYaml(YAML::Node node);
		Experiment::IAgent* getRandomIndividual();
		std::vector<Experiment::IAgent*> getMutatedChildren(
			const std::vector<Experiment::IAgent*>& parents);

	protected:
		bool m_initialized;
		CGP_TYPE m_agentType;

		// This is only used by standard CGP to indicate the type of values
		// the nodes process on. For now, that is limited to double and bool:
		std::string m_nodeProcessingType;

		YAML::Node m_mutationParameters;
		YAML::Node m_cgpGeneralParameters;

		// Each CGP type has different parameters specific to it. Leave those in the
		// node it came from for easy access:
		YAML::Node m_versionSpecificParameters;	

		// Get different, specific individual types:
		Experiment::IAgent* getRandomCGPFFANNIndividual();
		Experiment::IAgent* getRandomStandardCGPIndividual_double();
		Experiment::IAgent* getRandomStandardCGPIndividual_bool();

	private:
		virtual std::string classname() { return "CGPGenerator"; }
		
	};
}

