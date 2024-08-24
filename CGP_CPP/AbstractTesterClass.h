#pragma once

#include <vector>
#include "StandardCGPIndividual.h"

/*
* This class represents the interface that must be provided to the
* GeneralCGPSolver if a full class must be provided rather than a
* single testing function. It is meant to handle situations where a single
* testing function would be too limited or would waste significant time
* re-initializing the testing environment.
*/
namespace CGP
{
	class AbstractTesterClass
	{
	private:
		virtual std::string classname() { return "AbstractTesterClass"; }

	protected:
		AbstractTesterClass();
		virtual ~AbstractTesterClass();

	public:
		std::vector<double> virtual getPopulationScores(
			std::vector<CGP::AbstractCGPIndividual* > population,
			int fitnessTestRepeats, unsigned int maxStepsPerRun) = 0;
	};

}; // End CGP Namespace