#pragma once
#include <vector>
#include <unordered_map>
#include <tuple>
#include "Enums.h"
#include "IAgent.h"
// #include <json.hpp>  # https://github.com/nlohmann/json

namespace CGP
{
	/********************************************************************
	* Class: AbstractCGPIndividual
	* Purpose: This is the abstract individual from which all individuals
	*          must inherit.
	*/

		/********************************************************************
		* Name:
		* Parameters:
		* Returns:
		* Purpose:
		*********************************************************************/

	class AbstractCGPIndividual : public Experiment::IAgent
	{
	public:
		// Base class destructors must be virtual to allow derived classes to cleanup:
		virtual ~AbstractCGPIndividual();

	private:
		bool m_baseClassInitialized;

		virtual std::string classname() { return "AbstractCGPIndividual"; }

	protected:
		int m_inputSize;
		int m_outputSize;
		int m_rows;
		int m_columns;
		int m_columnsBack;
		int m_columnsForward;

		AbstractCGPIndividual();
		AbstractCGPIndividual(int inputSize, int outputSize, int rows, int columns,
			int colsBack, int colsForward);		

		/********************************************************************
		* Name: initializeAbstractIndividual
		* Parameters:
		*   - int inputSize: The number of inputs to expect.
		*   - int outputSize: The number of outputs to provide.
		*   - int rows: The number of rows of processing nodes.
		*   - int columns: The number of columns of processing nodes.
		*   - int colsBack: How many columns backward a node is allowed to reach for inputs.
		*   - int colsForward: How many columns forward a node is allowed to reach for inputs.
		* Returns: Void
		* Purpose: To avoid needing parameters to the constructor, needed
		*          variables will be passed to the initializer function.
		*********************************************************************/
		void initializeAbstractIndividual(int inputSize, int outputSize, int rows,
			int columns, int colsBack, int colsForward);

	public:
		/*****************************************************************
		* Name: randomize
		* Parameters: None
		* Returns: Void
		* Purpose: randomize modifies itself to be a random individual of its
		*          specified type. This isn't used as mutation, but rather to
		*          generate an individual with no parent.
		******************************************************************/
		virtual void randomize() = 0;

		/********************************************************************
		* Name: calculateOutputs
		* Parameters:
		*   - inputs: std::vector of int or double representing the problem or
		*             environmental inputs to be run through calculation.
		* Returns: The individual's calculated output as a std::vector of double
		* Purpose: Pass the inputs through the function the CGP individual
		*          has cobbled together and return its output.
		*          Child classes only need to provide the version that uses
		*          double values.
		*********************************************************************/
		// virtual std::vector<double> calculateOutputs(std::vector<double> inputs) = 0;
		// virtual std::vector<double> calculateOutputs(std::vector<int> inputs);

		/********************************************************************
		* Name: getOneMutatedChild
		* Parameters:
		*   - strategy: An enum representing the mutation strategy.
		*   - parameters: std::unordered_map of that strategy's parameters.
		* Returns: An AbstractCGPIndividual* which is a mutated version of this instance.
		* Purpose: Create a mutated individual using this instance as a parent.
		*********************************************************************/
		virtual AbstractCGPIndividual* getOneMutatedChild(
			MUTATION_STRATEGY strategy, std::unordered_map<std::string, double> parameters) = 0;


		/********************************************************************
		* Name: performOncePerEpochUpdates
		* Parameters:
		*   - population: std::vector of AbstractCGPIndividual* representing
		*                 the full population available in this current experiment.
		*   - epochScore: The score the best individual achieved during the last
		*                 epoch.
		* Returns: Void
		* Purpose: Make any updates that need to be once per epoch that
		*          can't be handled by each individual on their own.
		*********************************************************************/
		virtual void performOncePerEpochUpdates(
			std::vector<AbstractCGPIndividual*> population, double epochScore) = 0;

		/********************************************************************
		* Name: calculateOutputs
		* Parameters:
		*   - observation: std::vector of values. If the value type doesn't
		*                  make sense for the individual, it doesn't need to
		*                  implement that function.
		* Returns: A value of the same type as the vector input.
		* Purpose: Calculate the individual's output given an observation
		*********************************************************************/
		virtual std::vector<double> calculateOutputs(const std::vector<double>& observation);
		virtual std::vector<int> calculateOutputs(const std::vector<int>& observation);
		virtual std::vector<bool> calculateOutputs(const std::vector<bool>& observation);


		/********************************************************************
		* Name: printGenotype
		* Parameters: None
		* Returns: None
		* Purpose: Print a human-readable genotype to standard screen output.
		*********************************************************************/
		virtual void printGenotype() = 0;

		/********************************************************************
		* Name: constrain
		* Parameters:
		*   - value: The value output by any one of this individual's
		*            calculation nodes.
		* Returns: A double representing that value constrained to acceptable
		*          bounds.
		* Purpose: To constrain node outputs to acceptable values during calculations.
		*********************************************************************/
		virtual double constrain(double value);

		/********************************************************************
		* Name: resetForNewTimeSeries
		* Parameters: None
		* Returns: Void
		* Purpose: Perform any changes or updates that the class requires
		*          in order to ready itself for a new time series (if working
		*          in a reinforcement learning environment) or problem (if
		*          trying to solve a simple equation.
		*********************************************************************/
		virtual void resetForNewTimeSeries();

		/********************************************************************
		* Name: getPercentageNodesUsed
		* Parameters: None
		* Returns: A double between 0.0 and 100.0.
		* Purpose: For data collection, determine the percentage of the
		*          processing nodes that are being utilized by the current
		*          node layout.
		*********************************************************************/
		virtual double getPercentageNodesUsed() = 0;

		/********************************************************************
		* Name: saveSelfToJson
		* Parameters: 
		* Returns:
		* Purpose:
		*********************************************************************/
		// virtual void saveSelfToJson(std::vector<std::string> funcNames) = 0;
				
		/********************************************************************
		* Name: getColumnNumber
		* Parameters:
		*   - int nodeNumber: The node's place in the list.
		* Returns: int representing the node's column number.
		* Purpose: Calculate a node's column number given the layout of this
		*          individual.
		*********************************************************************/
		int getColumnNumber(int nodeNumber);

		/********************************************************************
		* Name: getValidInputNodeNumber
		* Parameters:
		*   - int nodeNumber: The number of the node which needs an input.
		* Returns: An integer of a valid node number.
		* Purpose: Given the node and how far forward and backward this CGP
		*          individual is allowed to search, randomly choose an input.
		*********************************************************************/
		int getValidInputNodeNumber(int nodeNumber);

		/********************************************************************
		* Name: getNodeNumberRange
		* Parameters:
		*   - int columnNumber: The column number for which we want the node
		                        number range.
		* Returns: An integer tuple of the min and max node numbers in that column.
		* Purpose: Determine the range of node numbers in a given column.
		*********************************************************************/
		std::tuple<int, int> getNodeNumberRange(int columnNumber);		
	};

}  // End CGP namespace
