#include "pch.h"
#include "GymWorldGenerator.h"
#include "Enums.h"
#include "GymTester.h"
#include "SimpleRepeat.h"
#include "DungeonRoomEnv.h"

#include "yaml-cpp/yaml.h"

#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>

extern unsigned int DEBUG_LEVEL;

namespace CGP
{
	GymWorldGenerator::GymWorldGenerator() :
		m_initialized(false),
		m_gymType(GYM_TYPE::UNDEFINED),
		m_environmentName(""),
		m_observationSize(-1),
		m_actionSize(-1),
		m_useArgMax(false),
		m_stepsPerRun(10000),
		m_runsPerIndividual(-1),
		m_runsPerIndividualToGrade(-1),
		m_fitnessCollapseFunction(FITNESS_COLLAPSE_FUNCTION::UNDEFINED),
		m_CGym1D(nullptr),
		m_CGym2D(nullptr)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;
	}

	GymWorldGenerator::~GymWorldGenerator()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		if (m_CGym1D != nullptr)
			delete m_CGym1D;

		if (m_CGym2D != nullptr)
			delete m_CGym2D;
	}

	bool GymWorldGenerator::loadConfigurationFromYaml(YAML::Node node)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		// Assume success until we find failure:
		bool configSuccess = true;
		
		// Everything must be under "InstanceWorld":
		if (node["InstanceWorld"])
		{
			m_config = node["InstanceWorld"];

			if (m_config["Type"].as<std::string>() == "PyGym")
			{
				m_gymType = GYM_TYPE::PyGYM;

			}
			else if (m_config["Type"].as<std::string>() == "CGym_1D")
			{
				m_gymType = GYM_TYPE::CGYM_1D;
			}
			else if (m_config["Type"].as<std::string>() == "CGym_2D")
			{
				m_gymType = GYM_TYPE::CGYM_2D;
			}
			else
			{
				std::cout << "Unrecognized Gym-Type: " << m_config["Type"].as<std::string>() << std::endl;
				configSuccess = false;
			}

			// If we've succeeded so far, carry on:
			if (configSuccess)
			{
				YAML::Node gymConfig;
				if (m_gymType == GYM_TYPE::PyGYM)
				{
					gymConfig = m_config["PyGym"];
				}
				else if (m_gymType == GYM_TYPE::CGYM_1D)
				{
					gymConfig = m_config["CGym_1D"];
				}
				else if (m_gymType == GYM_TYPE::CGYM_2D)
				{
					gymConfig = m_config["CGym_2D"];
				}

				// If we got a valid gym config, grab data from it:
				if (!gymConfig.IsNull())
				{
					m_observationSize = gymConfig["ObservationSize"].as<int>();
					m_actionSize = gymConfig["ActionSize"].as<int>();
					m_useArgMax = (gymConfig["ActionType"].as<std::string>() == "Discrete");
					m_stepsPerRun = static_cast<unsigned int>(gymConfig["MaxSteps"].as<int>());
					m_runsPerIndividual = m_config["IndividualGrading"]["NumberOfRuns"].as<int>();
					m_runsPerIndividualToGrade = m_config["IndividualGrading"]["RunsToGrade"].as<int>();
					m_environmentName = gymConfig["Name"].as<std::string>();
					
					// Get the fitness collapse function:
					std::string fitString = m_config["IndividualGrading"]["FitnessCollapseFunction"].as<std::string>();
					if (fitString == "Mean")
						m_fitnessCollapseFunction = FITNESS_COLLAPSE_FUNCTION::MEAN;
					else if (fitString == "Median")
						m_fitnessCollapseFunction = FITNESS_COLLAPSE_FUNCTION::MEDIAN;
					else if (fitString == "AvgOfMeanAndMedian")
						m_fitnessCollapseFunction = FITNESS_COLLAPSE_FUNCTION::AVG_OF_MEAN_AND_MEDIAN;
					else if (fitString == "MinOfMeanAndMedian")
						m_fitnessCollapseFunction = FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN;
					else if (fitString == "Minimum")
						m_fitnessCollapseFunction = FITNESS_COLLAPSE_FUNCTION::MINIMUM;
					else if (fitString == "Maximum")
						m_fitnessCollapseFunction = FITNESS_COLLAPSE_FUNCTION::MAXIMUM;
					else if (fitString == "MinOfMeanAndMedianPlusMin")
						m_fitnessCollapseFunction = FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN_PLUS_MIN;
					else
					{
						std::cout << "Unrecognized fitness collapse function: " << fitString << std::endl;
						configSuccess = false;
					}
				}
			}
		}
		else
		{
			configSuccess = false;
		}

		return configSuccess;		
	}

	bool GymWorldGenerator::createTestWorld()
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		bool success = true;
		if (m_gymType == GYM_TYPE::PyGYM)
		{
			GymTester::getInstance()->initialize(m_environmentName, m_fitnessCollapseFunction,
				m_useArgMax, m_observationSize, m_actionSize);
		}
		else if (m_gymType == GYM_TYPE::CGYM_1D)
		{
			// Check that we got one of the names we recognize:
			if (m_environmentName == "SimpleRepeat")
			{
				// The input/outpus size must be the same on a repeat-environment.
				// The number of times to test each individual for a single fitness value
				// must be provided in the configuration file:
				m_CGym1D = new CGym::SimpleRepeat(m_actionSize, m_config["CGym_1D"]["TestSteps"].as<int>());
			}
			else if (m_environmentName == "DungeonRoom")
			{
				// DungeonRoom has many parameters:
				YAML::Node params = m_config["CGym_1D"];
				
				m_CGym1D = new CGym::DungeonRoomEnv(
					params["MaxCoord"].as<int>(), // Size of the dungeon room
					params["PlayerAlwaysOnSouth"].as<bool>(), // If player always starts on the south wall.
					params["ScoreLoss_HitWall"].as<double>(),
					params["ScoreLoss_PerTurn"].as<double>(),
					params["ScoreLoss_DistanceFromExit"].as<double>(),
					params["MaxSteps"].as<int>(),
					params["KeyScore"].as<double>(), // Score gain for grabbing the key
					params["ExitScore"].as<double>(), // Score gain for reaching the exit
					params["RequireKeyToExit"].as<bool>(),
					params["MultipleRooms"].as<bool>(), // More than one room if exit is found
					params["DistanceToWallInObservation"].as<bool>(), // Adds 4 values to the obs-size
					params["Allow8DirectionMovement"].as<bool>(), // 4 or 8 direction movement.
					params["RepeatedKeys"].as<bool>() // If true, finding the key places another.
				);
			}
			else
			{
				std::cout << "Unrecognized CGym_1D name: " << m_environmentName << std::endl;
				success = false;
			}
		}
		else if (m_gymType == GYM_TYPE::CGYM_2D)
		{			
			// No 2D CGyms yet.
			std::cout << "Unrecognized CGym_2D name: " << m_environmentName << std::endl;
			success = false;
		}
		else
		{
			std::cout << "Unrecognized gym type." << std::endl;
			success = false;
		}	
		
		return success;
	}

	std::vector<double> GymWorldGenerator::getPopulationFitnesses(
		const std::vector<Experiment::IAgent*>& population)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		if (m_gymType == GYM_TYPE::PyGYM)
		{
			// Turn our vector of IAgent into a vector of AbstractCGPIndividual, as required by GymTester:
			std::vector<AbstractCGPIndividual*> testPop;

			for (unsigned int i = 0; i < population.size(); ++i)
				testPop.push_back(static_cast<AbstractCGPIndividual*>(population[i]));

			// Return the calculated fitnesses:
			// TODO: Implement a difference between runs and graded runs at this level and the
			//       GymTester level.
			return GymTester::getInstance()->getPopulationScores(testPop, m_runsPerIndividual, m_stepsPerRun);
		}
		else if (m_gymType == GYM_TYPE::CGYM_1D)
		{
			if (m_CGym1D != nullptr)
			{
				return getCGym1DPopulationFitnesses(population);
			}
			else
			{
				std::cout << "Unknown error. CGym_1D pointer is null." << std::endl;
				assert(false);
				return std::vector<double>{-1.0};  // Prevent warnings about not returning values.
			}
		}
		else if (m_gymType == GYM_TYPE::CGYM_2D)
		{
			if (m_CGym2D != nullptr)
			{
				return getCGym2DPopulationFitnesses(population);
			}
			else
			{
				std::cout << "Unknown error. CGym_2D pointer is null." << std::endl;
				assert(false);
				return std::vector<double>{-1.0};  // Prevent warnings about not returning values.
			}
		}
		else
		{
			std::cout << "Unknown GymType." << std::endl;
			assert(false);
			return std::vector<double>{-1.0}; // Prevent warnings about not returning values.
		}
	}

	std::vector<double> GymWorldGenerator::getCGym1DPopulationFitnesses(const std::vector<Experiment::IAgent*>& population)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		std::vector<double> fullPopulationFitnesses;
		std::vector<double> singleIndividualFitnesses;

		for (unsigned int i = 0; i < population.size(); ++i)
		{
			// Gather all fitnesses for a single individual:
			singleIndividualFitnesses.clear();
			for (int j = 0; j < m_runsPerIndividual; ++j)
			{
				singleIndividualFitnesses.push_back(getCGym1DSingleFitness(population[i]));
			}

			// Collapse those values down to a single number and add it to our population vector:
			fullPopulationFitnesses.push_back(applyFitnessCollapseFunction(singleIndividualFitnesses));
		}

		// Return the full population of fitness vectors:
		return fullPopulationFitnesses;	
	}

	double GymWorldGenerator::applyFitnessCollapseFunction(std::vector<double>& fitnesses)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		if (m_fitnessCollapseFunction == FITNESS_COLLAPSE_FUNCTION::UNDEFINED)
		{
			std::cout << "Undefined fitness collapse function." << std::endl;
			return -1.0;
		}

		// Single value? Return it:
		if (fitnesses.size() == 1)
			return fitnesses[0];
		
		std::sort(fitnesses.begin(), fitnesses.end());

		double retVal = -1.0;
		double min = fitnesses[0];
		double max = fitnesses[fitnesses.size() - 1];
		double mean = std::accumulate(fitnesses.begin(), fitnesses.end(), 0.0) / double(fitnesses.size());
		double median = fitnesses[fitnesses.size() / 2];
		
		// If the length is even, median calculation changes:
		if (fitnesses.size() % 2 == 0)
		{
			median = (median + fitnesses[(fitnesses.size() / 2) - 1]) / 2.0;
		}
		
		switch (m_fitnessCollapseFunction)
		{
			case FITNESS_COLLAPSE_FUNCTION::MEAN:
				retVal = mean;
			break;

			case FITNESS_COLLAPSE_FUNCTION::MEDIAN:
				retVal = median;
			break;

			case FITNESS_COLLAPSE_FUNCTION::AVG_OF_MEAN_AND_MEDIAN:
				retVal = (mean + median) / 2.0;
			break;

			case FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN:
				retVal = std::min(mean, median);
			break;

			case FITNESS_COLLAPSE_FUNCTION::MINIMUM:
				retVal = min;
			break;

			case FITNESS_COLLAPSE_FUNCTION::MAXIMUM:
				retVal = max;
			break;

			case FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN_PLUS_MIN:
				retVal = std::min(mean, median) + min;
			break;

			// Shouldn't happen:
			case FITNESS_COLLAPSE_FUNCTION::UNDEFINED:
				assert(false);
			break;
		}

		return retVal;
	}

	double GymWorldGenerator::getCGym1DSingleFitness(Experiment::IAgent* agent)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		bool done = false;
		double stepReward = 0.0;
		double fullReward = 0.0;
		std::vector<double> action;
		std::vector<double> observation;
		std::vector<std::string> info;

		observation = m_CGym1D->reset();

		while (!done)
		{
			action = agent->calculateOutputs(observation);

			m_CGym1D->step(action, observation, stepReward, done, info);
			fullReward += stepReward;
		}

		return fullReward;
	}

	std::vector<double> GymWorldGenerator::getCGym2DPopulationFitnesses(const std::vector<Experiment::IAgent*>& population)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		std::cout << "2D-CGym not yet implemented." << std::endl;
		return std::vector<double>{-1.0};
	}

	double GymWorldGenerator::getCGym2DSingleFitness(Experiment::IAgent* agent)
	{
		if (DEBUG_LEVEL > 4)
			std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		std::cout << "2D-CGym not yet implemented." << std::endl;
		return -1.0;
	}
	
}
