#pragma once
#include "AbstractCGPIndividual.h"
#include "BasicGene.h"
#include "CGPFunctions.h"
#include <functional>
#include <vector>
#include "JsonLib/json.hpp"
using json = nlohmann::json;

namespace CGP
{
	// The neuron that will be used by the CGPFFANNIndividual class:
	struct CGPNeuron
	{
		GENE_TYPE m_geneType;
		std::vector<int> m_inputs;
		std::vector<double> m_weights;
		int m_function;
		double m_preBias;
		double m_postBias;

		CGPNeuron(GENE_TYPE geneType)
		: m_geneType(geneType), m_function(-1), m_preBias(0.0), m_postBias(0.0) {}
	};

    class CGPFFANNIndividual :
        public AbstractCGPIndividual
    {
	public:
		CGPFFANNIndividual(int inputSize, int outputSize, int rows, int columns,
			int colsBack, double minWeight, double maxWeight, int minNeuronInputCount,
			int maxNeuronInputCount, bool useConstraint, double minConstraint,
			double maxConstraint, double minPreBias, double maxPreBias,
			double minPostBias, double maxPostBias,
			std::vector<std::string> functionList);

		virtual void randomize();
			
		virtual AbstractCGPIndividual* getOneMutatedChild(
			MUTATION_STRATEGY strategy, std::unordered_map<std::string, double> parameters);

		virtual void performOncePerEpochUpdates(
			std::vector<AbstractCGPIndividual*> population, double epochScore);
				
		virtual void printGenotype();				
				
		virtual void resetForNewTimeSeries();
				
		virtual double getPercentageNodesUsed();

		std::vector<double> calculateOutputs(const std::vector<double>& input);

		void writeSelfToJson(json& j);
		void readSelfFromJson(json& j);

	protected:
		double m_minWeight;
		double m_maxWeight;
		// Pre-bias is added before the NN function is applied.
		double m_minPreBias;
		double m_maxPreBias;
		bool m_mutatePreBias;
		// Post-bias is added after the NN function is applied:
		double m_minPostBias;
		double m_maxPostBias;
		bool m_mutatePostBias;
		bool m_useConstraint;
		double m_minConstraint;
		double m_maxConstraint;
		unsigned int m_minInputCount;
		unsigned int m_maxInputCount;
		std::vector<CGPNeuron> m_genotype;
		// Keeping our functions as a list of strings allows us to read/write
		// them to/from files:
		std::vector<std::string> m_functionStringList;
		// m_activeFunctions holds which functions can be selected during mutation.
		// This allows the available functions to be changed via mutation:
		std::vector<unsigned int> m_activeFunctions;

		std::vector<std::function<double(double)> > m_functionList;
		std::vector<int> m_activeGenes;
		
		virtual void calculateActiveGenes();
		
		virtual void mutateSelf(MUTATION_STRATEGY strategy,
			std::unordered_map<std::string, double> parameters);
		
		virtual void mutateSelf_probabilisticPerGene(
			std::unordered_map<std::string, double> parameters);
		
		virtual void mutateSelf_probabilisticPerValue(
			std::unordered_map<std::string, double> parameters);
		
		virtual void mutateSelf_activeGene(
			std::unordered_map<std::string, double> parameters);

		// If the difference between the minimum and maximum
		// bias is less than this value, that bias won't be
		// changed during mutation:
		static double BIAS_CHANGE_MINIMUM;

	private:
		virtual std::string classname() { return "CGPFFANNIndividual"; }
    };

}; // End CGP Namespace
