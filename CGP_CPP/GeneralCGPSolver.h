#pragma once
#include "AbstractCGPIndividual.h"
#include "AbstractTesterClass.h"
#include "Enums.h"
#include <functional>
#include <vector>
#include <mutex>
#include <queue>
#include <thread>
#include <chrono>
#include <unordered_map>

namespace CGP
{
	struct CGPSolverResults
	{
		int individualNumber;
		double fitness;

		CGPSolverResults(int ind, double fit)
			: individualNumber(ind), fitness(fit)
		{}
	};

	class GeneralCGPSolver
	{
	private:
		// An infinite-loop function designed to run in a thread and 
		// calculate individual fitnesses:
		void threadFitnessFunction();

		// Pass work to worker threads:
		std::mutex m_workQueueMutex;
		std::queue<int> m_workQueue;

		// Get results back from worker function:
		std::mutex m_resultsQueueMutex;
		std::queue<CGPSolverResults> m_resultsQueue;

		// threadFitnessFunction runs until it sees this value
		// as its input:
		static const int EXIT_THREAD = -1;
		
		// Our vector of threads:
		std::vector<std::thread*> m_threads;

		virtual std::string classname() { return "GeneralCGPSolver"; }

	protected:
		std::string getCurrentTimeString();

		// A function that will give us a random individual:		
		std::function<AbstractCGPIndividual* ()> m_randomIndividualGenerator;
		
		// A function to provide an individual's fitness. Higher 
		// values must always indicate a more fit individual:
		std::function<double(AbstractCGPIndividual*)> m_fitnessFunction;

		// If the fitness is difficult to calculate in a single function,
		// a class can be provided:
		AbstractTesterClass* m_fitnessClass;

		bool m_useFitnessClass;

		// Our population:
		std::vector<AbstractCGPIndividual*> m_population;
		
		// The maximum number of steps per run through the environment:
		unsigned int m_maxStepsPerRun;

		// Maximum number of epochs, even if we don't find a solution:
		unsigned long m_maxEpochs;

		// The number of individuals in each generation:
		unsigned int m_maxPopulationSize;

		// The number of threads to use when processing fitness:
		unsigned int m_processingThreadCount;

		// Output to the screen every X epochs:
		int m_screenOutputMod;

		// Score at which we stop processing and declare the problem solved:
		double m_solvedScore;

		// The number of times to test an individual to help mitigate the
		// impact of luck.
		unsigned int m_fitnessTestRepeats;

		// If True, the parent (best from last epoch) will be re-evaluated
		// in the current epoch. This will take more processing time, but
		// will force the parent to repeat their performance in order to stay
		// on top.  If environments with a lot of randomness, setting this value
		// to True can help prevent bad parents based on luck.
		bool m_reevaluateParent;

		// The function we'll use to turn a vector of fitnesses into a
		// single value that can be used for comparison:
		FITNESS_COLLAPSE_FUNCTION m_collapseFunction;

		// Saved to be queried after the experiment is over:
		AbstractCGPIndividual* m_finalParent;
		double m_finalParentScore;

		// How we mutate our individuals:
		MUTATION_STRATEGY m_mutationStrategy;
		std::unordered_map<std::string, double> m_mutationParameters;

		// Calculate all fitness values:
		void calculateAllFitness(std::vector<double>& fitnessValues);
		void calculateAllFitness_Function(std::vector<double>& fitnessValues);
		void calculateAllFitness_Class(std::vector<double>& fitnessValues);
		
		// This will sort the incoming vector as part of calculating fitness.
		// Don't use it if the order of the fitnesses matters:
		double applyFitnessCollapseFunction(std::vector<double>& fitVector);

	public:
		
		GeneralCGPSolver(						
			unsigned int maxStepsPerRun,
			unsigned long maxEpochs,
			unsigned int popSize,
			unsigned int threads,
			unsigned int fitnessRepeats,
			bool reevaluateParent,
			FITNESS_COLLAPSE_FUNCTION collapseFunction,
			int outputMod,
			double solvedScore,
			MUTATION_STRATEGY mutStrat,
		    std::unordered_map<std::string, double> mutParams,
			
			// Either a function to generate individuals or a starting
			// parent needs to be provided:
			std::function<AbstractCGPIndividual* ()> generator = nullptr,
			AbstractCGPIndividual* startingParent = nullptr,
			
			// Either a fitness function or a fitness class must be
			// provided:
			std::function<double(AbstractCGPIndividual*)> fitnessFunc = nullptr,
			AbstractTesterClass* fitnessClass = nullptr);
		~GeneralCGPSolver();

		// Run the full experiment:
		void runExperiment();

		// Get the final parent:
		AbstractCGPIndividual* getFinalParent() { return m_finalParent; }
		double getFinalScore() { return m_finalParentScore; }
	};
};
