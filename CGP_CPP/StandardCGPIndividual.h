#pragma once
#include "AbstractCGPIndividual.h"
#include "BasicGene.h"
#include <vector>
#include <functional>
#include "JsonLib/json.hpp"
using json = nlohmann::json;

namespace CGP
{
    template <class funcType> class StandardCGPIndividual :
        public AbstractCGPIndividual
    {
    protected:
        bool m_pRelevant;  // Some node/functions don't use P-processing.
        funcType m_minP;
        funcType m_maxP;
        funcType m_minConstraint;
        funcType m_maxConstraint;
        bool m_useConstraint;
        double m_minScale;
        double m_maxScale;
        bool m_useScale;
        FUNCTION_TYPE m_nodeFunctionType;
        // Keeping our functions as a list of strings allows us to read/write
        // them to/from files:
        std::vector<std::string> m_functionStringList;
        std::vector<std::function<funcType(funcType, funcType, funcType)> > m_functionList;        
        // m_activeFunctions holds which functions can be selected during mutation.
        // This allows the available functions to be changed via mutation:
        std::vector<unsigned int> m_activeFunctions;
        std::vector<BasicGene> m_genome;
        std::vector<int> m_activeGenes;

    public:                
        StandardCGPIndividual(int inputSize, int outputSize, int rows, int columns,
            int colsBack, bool pRelevant, funcType minP, funcType maxP, bool useConstraint,
            funcType minConstraint, funcType maxConstraint, bool useScale, 
            double minScale, double maxScale, FUNCTION_TYPE nodeFuncType,
            std::vector<std::string> inFunctionList);
        ~StandardCGPIndividual();

        // Defined by AbstractCGPIndividual:
		virtual void randomize();
		virtual std::vector<funcType> calculateOutputs(const std::vector<funcType>& inputs);
		virtual AbstractCGPIndividual* getOneMutatedChild(
			MUTATION_STRATEGY strategy, std::unordered_map<std::string, double> parameters);
		virtual void performOncePerEpochUpdates(
			std::vector<AbstractCGPIndividual*> population, double epochScore);
		virtual funcType constrain(funcType value);
        virtual void printGenotype();
		virtual void resetForNewTimeSeries();
        virtual double getPercentageNodesUsed();

        virtual void writeSelfToJson(json& j);
        virtual void readSelfFromJson(json& j);

        // Specific to StandardCGPIndividual:

        /********************************************************************
        * Name: calculateActiveGenes
        * Parameters: None
        * Returns: void
        * Purpose: Fill m_activeGenes with a list of all of all of the node
        *          numbers that are active in the current genome, in numerical
        *          order.
        *********************************************************************/
        virtual void calculateActiveGenes();
        virtual void mutateSelf(MUTATION_STRATEGY strategy,
                                std::unordered_map<std::string, double> parameters);
        virtual void mutateSelf_probabilisticPerGene(
            std::unordered_map<std::string, double> parameters);
        virtual void mutateSelf_probabilisticPerValue(
            std::unordered_map<std::string, double> parameters);
        virtual void mutateSelf_activeGene(
            std::unordered_map<std::string, double> parameters);
    private:
        virtual std::string classname() { return "StandardCGPIndividual"; }
    };
}

