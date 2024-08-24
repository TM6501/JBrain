#include "pch.h"
#include "GeneralCGPSolver.h"
#include "GymTester.h"
#include <mutex>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <chrono>  // Timestamps

const int CGP::GeneralCGPSolver::EXIT_THREAD;

extern unsigned int DEBUG_LEVEL;

namespace CGP
{
	GeneralCGPSolver::GeneralCGPSolver(		
		unsigned int maxStepsPerRun, unsigned long maxEpochs,
		unsigned int popSize, unsigned int threads,
		unsigned int fitnessRepeats, bool reevaluateParent,
		FITNESS_COLLAPSE_FUNCTION collapseFunction, int outputMod,
		double solvedScore, CGP::MUTATION_STRATEGY mutStrat,
		std::unordered_map<std::string, double> mutParams,
		std::function<AbstractCGPIndividual* ()> generator,
		AbstractCGPIndividual* startingParent,
		std::function<double(AbstractCGPIndividual*)> fitnessFunc,
		AbstractTesterClass* fitnessClass)
		: m_randomIndividualGenerator(generator),
		m_fitnessFunction(fitnessFunc),
		m_fitnessClass(fitnessClass),
		m_maxStepsPerRun(maxStepsPerRun),
		m_maxEpochs(maxEpochs),
		m_maxPopulationSize(popSize),
		m_processingThreadCount(threads),		
		m_screenOutputMod(outputMod),
		m_solvedScore(solvedScore),
		m_fitnessTestRepeats(fitnessRepeats),		
		m_reevaluateParent(reevaluateParent),
		m_collapseFunction(collapseFunction),
		m_finalParent(nullptr),
		m_finalParentScore(-1000.0),
		m_mutationStrategy(mutStrat),		 
		m_mutationParameters(mutParams)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Either the fitnessFunc or the FitnessClass must be defined:
		if (fitnessFunc == nullptr && fitnessClass == nullptr)
		{
			throw std::invalid_argument("Either a fitnessFunc or fitnessClass must be provided. Both cannot be null.");
		}

		// Need either a starting parent or a function to create individuals:
		if (generator == nullptr && startingParent == nullptr)
		{
			throw std::invalid_argument("Either a generator function or starting parent must be provided. Both cannot be null.");
		}

		if (startingParent != nullptr)
		{
			m_population.clear();
			m_population.push_back(startingParent);
		}

		// If a fitnessfunction is provided, we'll ignore the fitness class:
		m_useFitnessClass = (fitnessFunc == nullptr);

		// At least 1 processing thread is required:
		if (m_processingThreadCount < 1)
			m_processingThreadCount = 1;

		// Mod 0 causes issues. If they gave us 0, they probably meant they
		// didn't want any output:
		if (m_screenOutputMod == 0)
			m_screenOutputMod = -1;
	}

	GeneralCGPSolver::~GeneralCGPSolver()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;
	}

	std::string GeneralCGPSolver::getCurrentTimeString()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
		char buffer[50] = { 0 };
		std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%s", std::localtime(&now));

		return std::string(buffer);
	}

	void GeneralCGPSolver::threadFitnessFunction()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		int workValue;
		std::vector<double> allFitness;
		double fitness;
		bool foundWork;
		
		// Run forever:
		while (true)
		{
			foundWork = false;
						
			while (!foundWork)
			{
				// Get singular control of the work queue.
				m_workQueueMutex.lock();
				
				if (m_workQueue.empty())
				{
					m_workQueueMutex.unlock();
					
					// Sleep and wait for work:
					using namespace std::chrono_literals;
					std::this_thread::sleep_for(25ms);
				}
				else
				{
					workValue = m_workQueue.front();
					foundWork = true;

					// Leave the quit message in place if we found it:
					if (workValue != GeneralCGPSolver::EXIT_THREAD)
					{
						m_workQueue.pop();
						m_workQueueMutex.unlock();
					}
					else
					{
						m_workQueueMutex.unlock();
						return;  // Exit this function. We're done processing.
					}
				}
			} 

			allFitness.clear();
			for (unsigned int i = 0; i < m_fitnessTestRepeats; ++i)
				allFitness.push_back(m_fitnessFunction(m_population[workValue]));

			fitness = applyFitnessCollapseFunction(allFitness);

			// Create our return structure:
			CGPSolverResults x(workValue, fitness);

			// Put the results on the output queue:
			{
				std::lock_guard<std::mutex> outputLock(m_resultsQueueMutex);
				m_resultsQueue.push(x);
			}  // End braces will release m_resultsQueueMutex
		}
	}

	double GeneralCGPSolver::applyFitnessCollapseFunction(std::vector<double>& fitVector)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Only one value? Return it:
		if (fitVector.size() == 1)
		{
			return fitVector[0];
		}

		// For now, just go with min of median and mean:
		// Sort to get median:
		std::sort(fitVector.begin(), fitVector.end());

		double median = fitVector[fitVector.size() / 2];
		if (fitVector.size() % 2 == 0)
		{
			// If the length is even, need to include the value before:
			median = (median + fitVector[(fitVector.size() / 2) - 1]) / 2.0;
		}

		double mean = std::accumulate(fitVector.begin(), fitVector.end(), 0.0) / double(fitVector.size());

		// Return the minimum of the two:
		double retVal = median;
		if (mean < median)
			retVal = mean;

		return retVal;
	}


	void GeneralCGPSolver::calculateAllFitness(std::vector<double>& fitnessValues)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		if (m_useFitnessClass)
			calculateAllFitness_Class(fitnessValues);
		else
			calculateAllFitness_Function(fitnessValues);
	}

	void GeneralCGPSolver::calculateAllFitness_Class(std::vector<double>& fitnessValues)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		fitnessValues = m_fitnessClass->getPopulationScores(
		  m_population, m_fitnessTestRepeats, m_maxStepsPerRun);
	}

	void GeneralCGPSolver::calculateAllFitness_Function(std::vector<double>& fitnessValues)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// The calculation threads should be running. Put the work on their queues
		// and wait for it to be done:
		std::vector<bool> workComplete;
		fitnessValues.clear();

		// Record that no work has been done:
		for (unsigned int i = 0; i < m_population.size(); ++i)
		{
			workComplete.push_back(false);
			fitnessValues.push_back(-1.0);
		}

		// Put all of the requests on the work queue:		
		m_workQueueMutex.lock();
		for (unsigned int i = 0; i < m_population.size(); ++i)
			m_workQueue.push(i);
		m_workQueueMutex.unlock();

		// Gather results until all work is done:
		unsigned int resultsCollected = 0;
		while (resultsCollected < m_population.size())
		{
			m_resultsQueueMutex.lock();
			if (!m_resultsQueue.empty())
			{
				auto result = m_resultsQueue.front();
				m_resultsQueue.pop();
				m_resultsQueueMutex.unlock();

				fitnessValues[result.individualNumber] = result.fitness;
				workComplete[result.individualNumber] = true;
				++resultsCollected;
			}
			else
			{
				m_resultsQueueMutex.unlock();
				// Sleep and wait for work:
				using namespace std::chrono_literals;
				std::this_thread::sleep_for(25ms);
			}
		}
	}

	void GeneralCGPSolver::runExperiment()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// We don't want scientific notation on the output:
		std::cout.setf(std::ios_base::fixed);

		///////  ***********************
		// Use the starting population if it is there:
		/////// ************************

		// Fill our initial population:
		m_population.clear();
		for (unsigned int i = 0; i < m_maxPopulationSize; ++i)
			m_population.push_back(m_randomIndividualGenerator());

		// Create our fitness-processing threads, if we need them:
		m_threads.clear();
		if (!m_useFitnessClass)
		{
			for (unsigned int i = 0; i < m_processingThreadCount; ++i)
			{
				m_threads.push_back(
					new std::thread(&GeneralCGPSolver::threadFitnessFunction, this));
			}
		}

		std::vector<double> fitnesses;
		calculateAllFitness(fitnesses);
				
		// Find the max element:
		double parentScore = fitnesses[0];
		int maxIdx = 0;
		for (unsigned int i = 1; i < fitnesses.size(); ++i)
		{
			if (fitnesses[i] > parentScore)
			{
				parentScore = fitnesses[i];
				maxIdx = i;
			}
		}

		auto parent = m_population[maxIdx];

		for (unsigned int i = 0; i < m_maxEpochs; ++i)
		{
			if (m_screenOutputMod > 0 && i % m_screenOutputMod == 0)
			{
				// Output the current time and score:
				std::cout << getCurrentTimeString() 
					      << i
					      << " parent score: "
				          << std::setprecision(1) << parentScore 
					      << std::endl;
				/*std::cout << "  All scores: [";
				for (auto j = 0; j < fitnesses.size() - 1; ++j)
				{
					std::cout << std::setprecision(1) << fitnesses[j] << " ";
				}
				std::cout << std::setprecision(1) << 
					fitnesses[fitnesses.size() - 1] << "]" << std::endl;				*/
			}

			if (parentScore >= m_solvedScore)
			{
				if (m_screenOutputMod > 0)
				{
					std::cout << "Epoch " << i << " solution found." << std::endl;
					std::cout << "Final fitness: " << parentScore << std::endl;
				}
				break;
			}

			// Remove the parent from the population:
			m_population.erase(m_population.begin() + maxIdx);

			// Delete the rest of the individuals:
			for (unsigned int j = 0; j < m_population.size(); ++j)
				delete m_population[j];

			
			// Create the new population of individuals:
			m_population.clear();

			// Add the parent back into the population; force it to repeat its
			// performance if it wants to stay at the top of the pack.
			// This means that best-score WILL sometimes go down.
			if (m_reevaluateParent)
				m_population.push_back(parent);
						
			while (m_population.size() < m_maxPopulationSize)
			{
				m_population.push_back(
					parent->getOneMutatedChild(m_mutationStrategy, m_mutationParameters));
			}			

			// Get the next generation's fitnesses:
			calculateAllFitness(fitnesses);

			// Put the parent at the end if we aren't reevaluating it, so the
			// only way it is selected as the best is if it is better
			// (not just the same) as its progeny:
			if (!m_reevaluateParent)
			{
				m_population.push_back(parent);
				fitnesses.push_back(parentScore);
			}

			// Find the new parent and parent score:
			parentScore = fitnesses[0];
			maxIdx = 0;
			for (unsigned int j = 1; j < fitnesses.size(); ++j)
			{
				if (fitnesses[j] > parentScore)
				{
					parentScore = fitnesses[j];
					maxIdx = j;
				}
			}

			parent = m_population[maxIdx];
		}

		// Tell the processing threads to stop:
		m_workQueueMutex.lock();
		m_workQueue.push(GeneralCGPSolver::EXIT_THREAD);
		m_workQueueMutex.unlock();

		// Wait for all threads to finish:
		for (auto threadPointer : m_threads)
		{
			threadPointer->join();
			delete threadPointer;
		}
		m_threads.clear();

		// TODO: Store off any values the user of this class may later query //
		m_finalParent = parent;
		m_finalParentScore = parentScore;
			
		// Delete the individuals:
		for (auto indPointer : m_population)
		{
			// User of this class is now responsible for deleting the parent:
			if (indPointer != parent)
				delete indPointer;
		}
		m_population.clear();
	}
};
