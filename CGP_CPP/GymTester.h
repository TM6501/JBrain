#pragma once

#include "AbstractTesterClass.h"
#include <string>
#include <vector>
#include "Python.h"


#include <mutex>
#include <thread>
#include <queue>
#include "Enums.h"

/*
* This is a tester class designed to handle all the nonsense involved
* with using the Python-only OpenAI Gym to test the CGP individuals.
*/

namespace CGP
{	
	class GymTester : public AbstractTesterClass
	{
	private:
		// Constructor is private. Due to the nature of Python,
		// we have to ensure it is only called once:
		GymTester();

		virtual std::string classname() { return "GymTester"; }
				
	protected:
		// PyObjects we can create once and hold indefinitely:
		PyObject* m_gymModule;		
		PyObject* m_env;
		PyObject* m_envReset;
		PyObject* m_envStep;
		PyObject* m_envClose;
		PyObject* m_envRender;
		PyObject* m_gcModule;
		PyObject* m_gcCollectFunction;
		
		int m_observationSize;
		int m_actionSize;
		FITNESS_COLLAPSE_FUNCTION m_collapseFunction;
		bool m_useArgMax;
		bool m_initialized;
		std::string m_envName;

		void initializePython();
		void finalizePython();

		// double getIndScore(AbstractCGPIndividual* ind, int msRenderDelay = -1);
		double applyCollapseFunction(std::vector<double> indRewards);

	public:
		static GymTester* getInstance();

		void initialize(
			std::string envName,  // Name of the environment			
			FITNESS_COLLAPSE_FUNCTION fitFunc, // How we collapse many to a single value.
			bool useArgMax, // If true, give the env the idx of our max output
						   // rather than all of the outputs.
			int observationSize, // Size of observation the individual expects
			int actionSize // Size of action the environment will expect
		);

		~GymTester();

		std::vector<double> getPopulationScores(
			std::vector<CGP::AbstractCGPIndividual* > population,
			int runsPerFitness = 1, unsigned int maxStepsPerRun = 10000);

		// maxStopScore has been added to prevent potential
		// infinite loops related to an agent figuring out
		// how to keep an environment going forever.
		double getIndScore(AbstractCGPIndividual* ind, int msRenderDelay = -1,
			               int runsPerFitness = 1,
			               unsigned int maxStepsPerRun = 10000);

		// Render the individual's run to the screen and report back its
		// score on said run:
	    double renderIndividual(CGP::AbstractCGPIndividual* ind, int msRenderDelay = -1,
			                    unsigned int maxSteps = 10000);

	};
}; // End CGP namespace
