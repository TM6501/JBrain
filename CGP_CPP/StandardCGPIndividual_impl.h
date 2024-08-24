#pragma once
#include "pch.h"
#include "StandardCGPIndividual.h"
#include <random>
#include <typeinfo>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <numeric>

extern unsigned int DEBUG_LEVEL;

namespace CGP
{
    template <class funcType>
    StandardCGPIndividual<funcType>::StandardCGPIndividual(int inputSize,
        int outputSize, int rows, int columns, int colsBack,
        bool pRelevant, funcType minP, funcType maxP, bool useConstraint,
        funcType minConstraint, funcType maxConstraint, bool useScale,
        double minScale, double maxScale, FUNCTION_TYPE nodeFuncType,
        std::vector<std::string> inFunctionList)
        // 0 is passed for columnsForward. Currently, forward inputs
        // aren't supported since they can lead to circular dependencies:
        : AbstractCGPIndividual(inputSize, outputSize, rows, columns, colsBack, 0),
        m_pRelevant(pRelevant),
        m_minP(minP),
        m_maxP(maxP),
        m_minConstraint(minConstraint),
        m_maxConstraint(maxConstraint),
        m_useConstraint(useConstraint),
        m_minScale(minScale),
        m_maxScale(maxScale),
        m_useScale(useScale),
        m_nodeFunctionType(nodeFuncType),
        m_functionStringList(inFunctionList)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Reserve space for the genes. This doesn't create that many 
        // genes, it only means we won't reallocate the vector
        // when we're filling the space:
        m_genome.reserve(m_inputSize + m_outputSize + (m_rows * m_columns));

        // Fill the list of functions:
        for (unsigned int i = 0; i < m_functionStringList.size(); ++i)
        {
            if (m_nodeFunctionType == FUNCTION_TYPE::BOOL3_BOOL)
                m_functionList.push_back(
                    CGPFunctions::boolIn_boolOut::getFuncFromString(
                        m_functionStringList[i]));
            else if (m_nodeFunctionType == FUNCTION_TYPE::DOUBLE3_DOUBLE)
                m_functionList.push_back(
                    CGPFunctions::doubleIn_doubleOut::getFuncFromString(
                        m_functionStringList[i]));
            else
            {
                std::cout << "Unacceptable node function type for standard CGP." << std::endl;
                exit(-1);
            }
        }

        // All functions are active at the start:
        for (unsigned int i = 0; i < m_functionList.size(); ++i)
            m_activeFunctions.push_back(i);
    }

    template <class funcType>
    StandardCGPIndividual<funcType>::~StandardCGPIndividual<funcType>()
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::randomize()
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Clear our genome:
        m_genome.clear();

        // Add in all of the input genes:
        for (int i = 0; i < m_inputSize; ++i)
        {
            // Add an input gene; most gene values are meaningless
            // for an input gene:
            m_genome.push_back(
                BasicGene(GENE_TYPE::INPUT, -1, -1, -1, 1.0, 1.0));
        }

        // Create our random number generators:        
        static std::random_device rd;
        static std::mt19937 rng(rd());
        std::uniform_int_distribution<int> funcDistribution(0, int(m_activeFunctions.size()) - 1);
        std::uniform_real_distribution<double> pDistribution(m_minP, m_maxP);
        std::uniform_real_distribution<double> scaleDistribution(m_minScale, m_maxScale);

        // Add all of the processing genes:
        for (int i = 0; i < (m_rows * m_columns); ++i)
        {
            auto inX = getValidInputNodeNumber(i + m_inputSize);
            auto inY = getValidInputNodeNumber(i + m_inputSize);
            auto inFuncNum = funcDistribution(rng);
            // If P and Scale are used, they will always be real-valued:
            auto inP = pDistribution(rng);
            auto inScale = scaleDistribution(rng);

            m_genome.push_back(BasicGene(
                GENE_TYPE::PROCESSING,
                inX, inY, m_activeFunctions[inFuncNum], inP, inScale));
        }

        // Add all of the output genes. Their only concern is the X
        // input, from which they will pull their value:
        for (int i = 0; i < m_outputSize; ++i)
        {
            m_genome.push_back(BasicGene(
                GENE_TYPE::OUTPUT,
                getValidInputNodeNumber(i + m_inputSize + (m_rows * m_columns)),
                -1, -1, 1.0, 1.0));
        }

        // This should now be the size of our genome:
        assert(m_genome.size() == static_cast<unsigned int>(m_inputSize + m_outputSize + (m_rows * m_columns)));

        // Always recalculate our active genes when we do something that
        // could change them:
        calculateActiveGenes();
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::calculateActiveGenes()
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Sets don't have the [] operator, so we'll fake some of the 
        // set's functionality as we ensure values aren't double-added:
        m_activeGenes.clear();

        for (int i = 1; i <= m_outputSize; ++i)
        {
            m_activeGenes.push_back(int(m_genome.size() - i));
        }

        // Go through all genes already added and add those they require.
        // m_activeGenes will grow during this loop:
        std::vector<int> geneNumsToAdd = {};
        for (unsigned int i = 0; i < m_activeGenes.size(); ++i)
        {
            auto gene = m_genome[m_activeGenes[i]];
            geneNumsToAdd.clear();

            // Nothing to do with input genes:
            if (gene.type == GENE_TYPE::INPUT)
                continue;

            else if (gene.type == GENE_TYPE::OUTPUT)
            {
                geneNumsToAdd.push_back(gene.X);
            }

            else // Processing, add X and Y
            {
                geneNumsToAdd.push_back(gene.X);
                geneNumsToAdd.push_back(gene.Y);
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

    template <class funcType>
    std::vector<funcType> StandardCGPIndividual<funcType>::calculateOutputs(
        const std::vector<funcType>& inputs)
    {
        // if (DEBUG_LEVEL > 4)
        //    std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // From lowest to highest values, fill in the calculated values
        // for each node:
        std::unordered_map<int, funcType> nodeOutputs = {};

        for (unsigned int i = 0; i < inputs.size(); ++i)
        {
            nodeOutputs.emplace(i, inputs[i]);
        }

        // Go through the rest of the active nodes. Generally, there is no need
        // to run calculations for most of the nodes:
        for (unsigned int i = 0; i < m_activeGenes.size(); ++i)
        {
            auto gene = m_genome[m_activeGenes[i]];

            // Ignore input genes:
            if (gene.type == GENE_TYPE::INPUT)
                continue;

            // Just grab the X value if it is an output gene:
            else if (gene.type == GENE_TYPE::OUTPUT)
            {
                // The value should be there:
                assert(nodeOutputs.find(gene.X) != nodeOutputs.end());

                nodeOutputs.emplace(m_activeGenes[i], nodeOutputs[gene.X]);
            }

            // Do processing:
            else
            {
                // The input values should be there:
                assert(nodeOutputs.find(gene.X) != nodeOutputs.end());
                assert(nodeOutputs.find(gene.Y) != nodeOutputs.end());

                auto inX = nodeOutputs[gene.X];
                auto inY = nodeOutputs[gene.Y];

                funcType calculatedOutput;
                if (m_pRelevant)
                    calculatedOutput = static_cast<funcType>(
                        m_functionList[gene.F](inX, inY,
                            static_cast<funcType>(gene.P)));
                else
                    calculatedOutput = static_cast<funcType>(
                        m_functionList[gene.F](inX, inY, inY));

                // Apply scaling if needed:
                if (m_useScale && m_nodeFunctionType != FUNCTION_TYPE::BOOL3_BOOL)
                    calculatedOutput = static_cast<funcType>(calculatedOutput * gene.scale);

                // Constrain the value:
                calculatedOutput = constrain(calculatedOutput);

                // Put the value into the node outputs:
                nodeOutputs.emplace(m_activeGenes[i], calculatedOutput);
            }
        }

        // Push all of the output values to our return vector:
        std::vector<funcType> retValue;

        for (unsigned int i = static_cast<int>(m_genome.size()) - m_outputSize;
            i < m_genome.size(); ++i)
        {
            retValue.push_back(nodeOutputs[i]);
        }

        return retValue;
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::mutateSelf(MUTATION_STRATEGY strategy,
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
            std::cout << "Non-implemented mutation strategy." << std::endl;
            assert(false);
        }

        calculateActiveGenes();
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::mutateSelf_probabilisticPerValue(
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
        std::uniform_int_distribution<int> funcDistribution(0, int(m_activeFunctions.size()) - 1);
        std::uniform_real_distribution<double> pDistribution(m_minP, m_maxP);
        std::uniform_real_distribution<double> scaleDistribution(m_minScale, m_maxScale);
        std::uniform_real_distribution<double> mutateChanceDistribution(0.0, 100.0);

        // Step through every value of every gene, mutating if the
        // random chance dictates that we should. Start after
        // the input genes that can't be mutated:
        for (unsigned int i = static_cast<unsigned int>(m_inputSize);
             i < m_genome.size(); ++i)
        {
            // Output gene; X is its only mutable value:
            if (m_genome[i].type == GENE_TYPE::OUTPUT)
            {
                if (mutateChanceDistribution(rng) < percentage)
                    m_genome[i].X = getValidInputNodeNumber(i);
            }

            // Processing node. Check each value independently:
            else
            {
                // X:
                if (mutateChanceDistribution(rng) < percentage)
                {
                    // Keep selecting new inputs until it changes.
                    // Exclude the extremely rare case of having only a single
                    // input and being the first processing node. This
                    // will create an infinite loop:
                    if (!(m_inputSize < 2 && i == 1))
                    {
                        int startX = m_genome[i].X;
                        while (startX == m_genome[i].X)
                            m_genome[i].X = getValidInputNodeNumber(i);
                    }

                }

                // Y:
                if (mutateChanceDistribution(rng) < percentage)
                {
                    if (!(m_inputSize < 2 && i == 1))
                    {
                        int startY = m_genome[i].Y;
                        while (startY == m_genome[i].Y)
                            m_genome[i].Y = getValidInputNodeNumber(i);
                    }
                }

                // F: Function change if there's at least 2 functions:
                if (m_activeFunctions.size() > 1 &&
                    mutateChanceDistribution(rng) < percentage)
                {
                    int startFunc = m_genome[i].F;

                    // Keep selecting new values until a new one appears:
                    while (startFunc == m_genome[i].F)
                        m_genome[i].F = m_activeFunctions[funcDistribution(rng)];
                }

                // No need to confirm P and Scale change; they are real-values:
                // P:
                if (m_pRelevant && mutateChanceDistribution(rng) < percentage)
                    m_genome[i].P = pDistribution(rng);

                // Scale:
                if (m_useScale && mutateChanceDistribution(rng) < percentage)
                    m_genome[i].scale = scaleDistribution(rng);
            }
        }
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::mutateSelf_probabilisticPerGene(
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
        std::uniform_int_distribution<int> funcDistribution(0, int(m_activeFunctions.size()) - 1);
        std::uniform_real_distribution<double> pDistribution(m_minP, m_maxP);
        std::uniform_real_distribution<double> scaleDistribution(m_minScale, m_maxScale);
        std::uniform_real_distribution<double> mutateChanceDistribution(0.0, 100.0);

        // Build the list of gene parts we'll mutate:
        std::vector<std::string> genePartList = { "X", "Y" };

        if (m_activeFunctions.size() > 1)
            genePartList.push_back("F");

        if (m_pRelevant)
            genePartList.push_back("P");

        if (m_useConstraint)
            genePartList.push_back("scale");

        std::uniform_int_distribution<int> genePartDistribution(0,
            static_cast<int>(genePartList.size()) - 1);

        // Step through every gene, mutate one of its values if the random
        // chance determines we should mutate:
        for (unsigned int i = static_cast<unsigned int>(m_inputSize);
             i < m_genome.size(); ++i)
        {
            // Check if we're going to mutate this gene:
            if (mutateChanceDistribution(rng) < percentage)
            {
                // Yes, we're mutating this gene. Determine what part:
                int genePart = genePartDistribution(rng);

                // If we're mutating an output gene, it has to be X.
                // Ignore the random selection:
                if (m_genome[i].type == GENE_TYPE::OUTPUT)
                    genePart = 0;

                // Mutate X:
                if (genePartList[genePart] == "X")
                {
                    // Keep selecting new inputs until it changes.
                    // Exclude the extremely rare case of having only a single
                    // input and being the first processing node. This
                    // will create an infinite loop:
                    if (!(m_inputSize < 2 && i == 1))
                    {
                        int startX = m_genome[i].X;
                        while (startX == m_genome[i].X)
                            m_genome[i].X = getValidInputNodeNumber(i);
                    }
                }

                // Mutate Y:
                else if (genePartList[genePart] == "Y")
                {
                    if (!(m_inputSize < 2 && i == 1))
                    {
                        int startY = m_genome[i].Y;
                        while (startY == m_genome[i].Y)
                            m_genome[i].Y = getValidInputNodeNumber(i);
                    }
                }

                // Mutate F:
                else if (genePartList[genePart] == "F")
                {
                    int startFunc = m_genome[i].F;

                    // Keep selecting new values until a new one appears:
                    while (startFunc == m_genome[i].F)
                        m_genome[i].F = m_activeFunctions[funcDistribution(rng)];
                }

                // Mutate P:
                else if (genePartList[genePart] == "P")
                {
                    m_genome[i].P = pDistribution(rng);
                }

                // Scale:
                else
                {
                    m_genome[i].scale = scaleDistribution(rng);
                }
            }
        }
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::mutateSelf_activeGene(
        std::unordered_map<std::string, double> parameters)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Make sure the values we expect are there:
        assert(parameters.find("geneCount") != parameters.end());

        // Get the percentage of values that we will modify:
        // Convert the double to an integer:
        int genesToModify = static_cast<int>(parameters.find("geneCount")->second + 0.49);
        int genesModified = 0;

        // Create our random number generators:
        static std::random_device rd;
        static std::mt19937 rng(rd());
        std::uniform_int_distribution<int> funcDistribution(0, int(m_activeFunctions.size()) - 1);
        std::uniform_real_distribution<double> pDistribution(m_minP, m_maxP);
        std::uniform_real_distribution<double> scaleDistribution(m_minScale, m_maxScale);
        std::uniform_int_distribution<int> geneSelectionDistribution(
            m_inputSize, m_inputSize + m_outputSize + (m_rows * m_columns) - 1);

        // Build the list of gene parts we'll mutate:        
        std::vector<std::string> genePartList = { "X", "Y" };

        if (m_activeFunctions.size() > 1)
            genePartList.push_back("F");

        if (m_pRelevant)
            genePartList.push_back("P");

        if (m_useConstraint)
            genePartList.push_back("scale");

        std::uniform_int_distribution<int> genePartDistribution(0,
            static_cast<int>(genePartList.size()) - 1);

        while (genesModified < genesToModify)
        {
            // Select the gene to modify:
            int toModify = geneSelectionDistribution(rng);

            // Determine the part of the gene to modify:
            // Yes, we're mutating this gene. Determine what part:
            int genePart = genePartDistribution(rng);

            // If we're mutating an output gene, it has to be X.
            // Ignore the random selection:
            if (m_genome[toModify].type == GENE_TYPE::OUTPUT)
                genePart = 0;

            // Mutate X:
            if (genePartList[genePart] == "X")
            {
                // Keep selecting new inputs until it changes.
                // Exclude the extremely rare case of having only a single
                // input and being the first processing node. This
                // will create an infinite loop:
                if (!(m_inputSize < 2 && toModify == 1))
                {
                    int startX = m_genome[toModify].X;
                    while (startX == m_genome[toModify].X)
                        m_genome[toModify].X = getValidInputNodeNumber(toModify);
                }
            }

            // Mutate Y:
            else if (genePartList[genePart] == "Y")
            {
                if (!(m_inputSize < 2 && toModify == 1))
                {
                    int startY = m_genome[toModify].Y;
                    while (startY == m_genome[toModify].Y)
                        m_genome[toModify].Y = getValidInputNodeNumber(toModify);
                }
            }

            // Mutate F:
            else if (genePartList[genePart] == "F")
            {
                int startFunc = m_genome[toModify].F;

                // Keep selecting new values until a new one appears:
                while (startFunc == m_genome[toModify].F)
                    m_genome[toModify].F = m_activeFunctions[funcDistribution(rng)];
            }

            // Mutate P:
            else if (genePartList[genePart] == "P")
            {
                m_genome[toModify].P = pDistribution(rng);
            }

            // Scale:
            else
            {
                m_genome[toModify].scale = scaleDistribution(rng);
            }

            // If this was an active gene, increment how many active genes
            // we've modified:
            if (std::find(m_activeGenes.begin(), m_activeGenes.end(), toModify)
                != m_activeGenes.end())
            {
                ++genesModified;
            }
        }
    }

    template <class funcType>
    AbstractCGPIndividual* StandardCGPIndividual<funcType>::getOneMutatedChild(
        MUTATION_STRATEGY strategy, std::unordered_map<std::string, double> parameters)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Create a new individual with our constructor arguments:
        AbstractCGPIndividual* retVal = new StandardCGPIndividual<funcType>(
            m_inputSize, m_outputSize, m_rows, m_columns, m_columnsBack,
            m_pRelevant, m_minP, m_maxP, m_useConstraint, m_minConstraint,
            m_maxConstraint, m_useScale, m_minScale, m_maxScale,
            m_nodeFunctionType, m_functionStringList);

        // Give that individual a copy of our genome and active functions.
        // The = operator on vectors makes a copy:
        static_cast<StandardCGPIndividual*>(retVal)->m_genome = m_genome;
        static_cast<StandardCGPIndividual*>(retVal)->m_activeGenes = m_activeGenes;
        static_cast<StandardCGPIndividual*>(retVal)->m_activeFunctions = m_activeFunctions;

        // Let the new individual mutate itself:
        static_cast<StandardCGPIndividual*>(retVal)->mutateSelf(
            strategy, parameters);

        return retVal;
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::performOncePerEpochUpdates(
        std::vector<AbstractCGPIndividual*> population, double epochScore)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;
    }

    template <class funcType>
    funcType StandardCGPIndividual<funcType>::constrain(funcType value)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Assume no constraints for now:
        funcType retValue = value;

        // Constrain if needed. Can't use fmin/fmax in case of different
        // value types:
        if (m_useConstraint)
        {
            if (retValue > m_maxConstraint)
                retValue = m_maxConstraint;

            if (retValue < m_minConstraint)
                retValue = m_minConstraint;
        }

        return retValue;
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::writeSelfToJson(json& j)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Standard to all CGP:
        j["individualType"] = "StandardCGPIndividual";
        j["inputSize"] = m_inputSize;
        j["outputSize"] = m_outputSize;
        j["genoTypeSize"] = m_genome.size();
        
        // P-related:
        j["pRelevant"] = m_pRelevant;
        if (m_pRelevant)
        {
            j["minP"] = m_minP;
            j["maxP"] = m_maxP;
        }
        
        // Constraints:
        if (m_nodeFunctionType != FUNCTION_TYPE::BOOL3_BOOL)
        {
            j["minConstraint"] = m_minConstraint;
            j["maxConstraint"] = m_maxConstraint;
        }

        // Scaling:
        j["useScale"] = m_useScale;
        j["minScale"] = m_minScale;
        j["maxScale"] = m_maxScale;

        // The type of node functions we use:
        j["nodeFunctionType"] = FunctionTypeToString(m_nodeFunctionType);

        // The list of functions as strings; this includes those that may
        // not be currently active:
        j["functionList"] = json::array();
        for (unsigned int i = 0; i < m_functionStringList.size(); ++i)
            j["functionList"][i] = m_functionStringList[i];

        // Active functions are those that can be selected by a mutation.
        // Inactive functions may still be in the genome.
        j["activeFunctions"] = json::array();
        for (unsigned int i = 0; i < m_activeFunctions.size(); ++i)
            j["activeFunctions"][i] = m_activeFunctions[i];

        // Our genome:
        j["genes"] = json::array();
        for (unsigned int i = 0; i < m_genome.size(); ++i)
        {
            // Need the type for every gene:
            j["genes"][i]["type"] = GeneTypeToString(m_genome[i].type);
            
            // All non-input need the X value:
            if (m_genome[i].type != GENE_TYPE::INPUT)
                j["genes"][i]["X"] = m_genome[i].X;
            
            // Processing is the only type that needs the rest:
            if (m_genome[i].type == GENE_TYPE::PROCESSING)
            {
                j["genes"][i]["Y"] = m_genome[i].Y;
                j["genes"][i]["F"] = m_genome[i].F;  // Index into function list

                if (m_pRelevant)
                    j["genes"][i]["P"] = m_genome[i].P;

                if (m_useScale)
                    j["genes"][i]["scale"] = m_genome[i].scale;
            }
        }
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::readSelfFromJson(json& j)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        std::string cgpType = j["individualType"];
        if (cgpType != "StandardCGPIndividual")
        {
            std::cout << "Wrong type! Can't read in " << cgpType << std::endl;
            return;
        }

        // Cleare the genome and allocate the size needed to avoid later,
        // costly re-allocations:
        int genotypeSize = j["genoTypeSize"].get<int>();
        m_genome.clear();
        m_genome.reserve(genotypeSize);

        m_inputSize = j["inputSize"].get<int>();
        m_outputSize = j["outputSize"].get<int>();

        // P-related:
        m_pRelevant = j["pRelevant"].get<bool>();
        if (m_pRelevant)
        {
            m_minP = j["minP"].get<double>();
            m_maxP = j["maxP"].get<double>();
        }

        // The type of node functions we use:
        m_nodeFunctionType = StringToFunctionType(j["nodeFunctionType"]);

        // Constraints don't matter if using booleans:
        if (m_nodeFunctionType != FUNCTION_TYPE::BOOL3_BOOL)
        {
            m_minConstraint = j["minConstraint"].get<double>();
            m_maxConstraint = j["maxConstraint"].get<double>();
        }

        // Scaling:
        m_useScale = j["useScale"].get<bool>();
        m_minScale = j["minScale"].get<double>();
        m_maxScale = j["maxScale"].get<double>();
        
        // The list of functions as strings; this includes those that may
        // not be currently active:
        m_functionStringList.clear();        
        for (unsigned int i = 0; i < j["functionList"].size(); ++i)
            m_functionStringList.push_back(j["functionList"][i].get<std::string>());            

        // Fill the list of functions:
        m_functionList.clear();
        m_functionList.reserve(m_functionStringList.size());
        for (unsigned int i = 0; i < m_functionStringList.size(); ++i)
        {
            if (m_nodeFunctionType == FUNCTION_TYPE::BOOL3_BOOL)
                m_functionList.push_back(
                    CGPFunctions::boolIn_boolOut::getFuncFromString(
                        m_functionStringList[i]));
            else if (m_nodeFunctionType == FUNCTION_TYPE::DOUBLE3_DOUBLE)
                m_functionList.push_back(
                    CGPFunctions::doubleIn_doubleOut::getFuncFromString(
                        m_functionStringList[i]));
            else
            {
                std::cout << "Unacceptable node function type for standard CGP." << std::endl;
                exit(-1);
            }
        }

        // Active functions are those that can be selected by a mutation.
        // Inactive functions may still be in the genome.
        m_activeFunctions.clear();
        for (unsigned int i = 0; i < j["activeFunctions"].size(); ++i)
            m_activeFunctions[i] = j["activeFunctions"][i].get<unsigned int>();

        // Determine the active genes for ourselves:
        calculateActiveGenes();
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::resetForNewTimeSeries()
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

    }

    template <class funcType>
    double StandardCGPIndividual<funcType>::getPercentageNodesUsed()
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        return (double(m_activeGenes.size()) / double(m_genome.size())) * 100.0;
    }

    template <class funcType>
    void StandardCGPIndividual<funcType>::printGenotype()
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

                    if (m_pRelevant)
                        std::cout << "   P:" << std::setprecision(3) << gene.P;

                    if (m_useScale)
                        std::cout << "   Scale:" << std::setprecision(3) << gene.scale;
                    
                    std::cout << "   F:" << m_functionStringList[gene.F];
                    std::cout << std::endl;
                }
            }
        }

        std::cout << "Active genes: ";
        for (int i : m_activeGenes)
            std::cout << i << " ";
    }
}
