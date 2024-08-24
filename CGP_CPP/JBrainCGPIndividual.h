#pragma once

/*
* This is the variation on the AbstractCGPIndividual that will be used by the
* JBrains for updating all variables within a neuron.
*/

#include "AbstractCGPIndividual.h"
#include "yaml-cpp/yaml.h"
#include <vector>
#include <string>
#include <functional>
#include "BasicGene.h"

using json = nlohmann::json;

namespace CGP
{
	class JBrainCGPIndividual : public AbstractCGPIndividual
	{
	public:
		JBrainCGPIndividual(
			unsigned int inputSize,
			unsigned int outputSize,
			unsigned int rows,
			unsigned int columns,
			unsigned int colsBack,
			float minP, float maxP,
			std::vector<std::string> functionStringList,
			std::vector<std::function<double(double, double, double)> > functionList,
			float minConstraint, float maxConstraint, bool useConstraint
		);
		~JBrainCGPIndividual();
		JBrainCGPIndividual(const JBrainCGPIndividual& other);

		void randomize() override;

		// Mutation is a function we must override, so it returns a base class
		// pointer allocated as a derived class instance:
		AbstractCGPIndividual* getOneMutatedChild(MUTATION_STRATEGY strategy,
	      std::unordered_map<std::string, double> parameters) override;
		
		void mutateSelf(MUTATION_STRATEGY strategy);

		// Once per epoch updates aren't required for brain individuals,
		// but the function must still exist:
		void performOncePerEpochUpdates(
			std::vector<AbstractCGPIndividual*> population, double epochScore) override;

		// Only the double-version needs to be overriden:
		std::vector<double> calculateOutputs(const std::vector<double>& observation) override;

		void printGenotype() override;
		void resetForNewTimeSeries() override;
		double getPercentageNodesUsed() override;

		// From the IAgent interface:
		void writeSelfToJson(json& j) override;
		void readSelfFromJson(json& j) override;

		// Operators to make copying/comparing easier:
		JBrainCGPIndividual& operator=(const JBrainCGPIndividual& rhs);
		bool operator==(const JBrainCGPIndividual& rhs);

		// Allocate a new JBrainCGPIndividual from JSON:
		static JBrainCGPIndividual* getCGPIndividualFromJson(json& j);

		void setMinConstraint(const float& val) { m_minConstraint = val; }
		void setMaxConstraint(const float& val) { m_maxConstraint = val; }
		void setMinP(const float& val) { m_minP = val; }
		void setMaxP(const float& val) { m_maxP = val; }

	protected:
		// All brain functions will be double-in-double-out:
		std::vector<std::string> m_functionStringList;
		std::vector<std::function<double(double, double, double)> > m_functionList;
		std::vector<NoScaleGene> m_genome;
		std::vector<int> m_activeGenes;

		// Determine which genes need to be calculated and which can be ignored
		// because their output is not used:
		void calculateActiveGenes();

		

	private:
		double constrain(double value) override;

		virtual std::string classname() { return "JBrainCGPIndividual"; }

		float m_minP;
		float m_maxP;
		float m_minConstraint;
		float m_maxConstraint;
		bool m_useConstraint;
	};

}  // End CGP namespace