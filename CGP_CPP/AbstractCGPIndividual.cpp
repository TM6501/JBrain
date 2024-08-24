#include "pch.h"
#include "AbstractCGPIndividual.h"
#include <assert.h>
#include <random>
#include <algorithm>  // min/max

// For debugging:
#include <iostream>

// extern unsigned int DEBUG_LEVEL;
// Not using here; too many low-level utility functions.

namespace CGP
{
	AbstractCGPIndividual::AbstractCGPIndividual()
		: m_baseClassInitialized(false), m_inputSize(-1), m_outputSize(-1), m_rows(-1),
		m_columns(-1), m_columnsBack(-1), m_columnsForward(-1)
	{
	}

	AbstractCGPIndividual::AbstractCGPIndividual(int inputSize, int outputSize,
		int rows, int columns, int colsBack, int colsForward)
	{
		initializeAbstractIndividual(inputSize, outputSize, rows, columns,
			colsBack, colsForward);
	}

	AbstractCGPIndividual::~AbstractCGPIndividual()
	{
	}

	void AbstractCGPIndividual::initializeAbstractIndividual(int inputSize, int outputSize,
		int rows, int columns, int colsBack, int colsForward)
	{
		m_baseClassInitialized = true;
		m_inputSize = inputSize;
		m_outputSize = outputSize;
		m_rows = rows;
		m_columns = columns;
		m_columnsBack = colsBack;
		m_columnsForward = colsForward;
	}

	//std::vector<double> AbstractCGPIndividual::calculateOutputs(std::vector<int> inputs)
	//{
	//	// Convert the vector of int to double, then call the double version
	//	std::vector<double> x = { 1.0 };
	//	return x;
	//}

	double AbstractCGPIndividual::constrain(double value)
	{
		// Standard behavior just returns the value. This function can
		// be overriden if different behavior is required.
		return value;
	}

	void AbstractCGPIndividual::resetForNewTimeSeries()
	{
		// Do nothing. Individuals can override if they have processing to do
		// to reset themselves.
	}

	int AbstractCGPIndividual::getColumnNumber(int nodeNumber)
	{
		// This function cannot be called until the base class is initialized:
		assert(m_baseClassInitialized);

		int totalLength = m_inputSize + m_outputSize + (m_rows * m_columns);

		int retVal = -1;

		// Input nodes are in column 0:
		if (nodeNumber < m_inputSize)
			retVal = 0;

		// Output nodes are in the final column:
		else if (nodeNumber >= totalLength - m_outputSize)
			retVal = m_columns + 1;

		// Otherwise, it is a processing node:
		else
			retVal = ((nodeNumber - m_inputSize) / m_rows) + 1;

		return retVal;
	}

	int AbstractCGPIndividual::getValidInputNodeNumber(int nodeNumber)
	{
		// This function cannot be called until the base class is initialized:
		assert(m_baseClassInitialized);

		// Get this node's column number:
		int colNum = getColumnNumber(nodeNumber);

		// Get the range of valid column numbers:
		int maxCol = std::min(colNum + m_columnsForward, m_columns);
		int minCol = std::max(colNum - m_columnsBack, 0);

		// Get the two ranges of nodes:
		int min1, min2, max1, max2;
		std::tie(min1, min2) = getNodeNumberRange(minCol);
		std::tie(max1, max2) = getNodeNumberRange(maxCol);

		// Choose a new random number, making sure we don't select ourselves:
		static std::random_device rd;
		static std::mt19937 rng(rd());
		std::uniform_int_distribution<int> distribution(min1, max2);

		int retVal = nodeNumber;
		while (retVal == nodeNumber)
			retVal = distribution(rng);

		return retVal;
	}

	std::tuple<int, int> AbstractCGPIndividual::getNodeNumberRange(int columnNumber)
	{
		// This function cannot be called until the base class is initialized:
		assert(m_baseClassInitialized);

		// Input node?
		if (columnNumber == 0)
			return std::make_tuple(0, m_inputSize - 1);

		// Output node?
		else if (columnNumber > m_columns)
		{
			int maxNode = m_inputSize + m_outputSize + (m_rows * m_columns) - 1;
			int minNode = maxNode - m_outputSize + 1;
			return std::make_tuple(minNode, maxNode);
		}

		// Processing node
		else
		{
			int minNode = ((columnNumber - 1) * m_rows) + m_inputSize;
			return std::make_tuple(minNode, minNode + (m_rows - 1));
		}			
	}

	std::vector<double> AbstractCGPIndividual::calculateOutputs(const std::vector<double>& observation)
	{
		std::cout << "In the abstract double calculate outputs." << std::endl;
		std::cout << "If this is to be used, the individual needs to implement it." << std::endl;
		return observation;
	}

	std::vector<int> AbstractCGPIndividual::calculateOutputs(const std::vector<int>& observation)
	{
		std::cout << "In the abstract int calculate outputs." << std::endl;
		std::cout << "If this is to be used, the individual needs to implement it." << std::endl;
		return observation;
	}

	std::vector<bool> AbstractCGPIndividual::calculateOutputs(const std::vector<bool>& observation)
	{
		std::cout << "In the abstract bool calculate outputs." << std::endl;
		std::cout << "If this is to be used, the individual needs to implement it." << std::endl;
		return observation;
	}

}  // End CGP namespace.

