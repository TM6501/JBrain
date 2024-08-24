#pragma once
#include "IInstanceWorld.h"
#include "IAgent.h"

#include "Enums.h"
#include "AbstractGymEnv.h"
#include "GymTester.h"

#include "yaml-cpp/yaml.h"

#include <string>
#include <vector>

namespace CGP {

    // The GymWorldGenerator will act as the interface between the
    // Experiment and the Gym instances. It will handle both the Python
    // and CGym Gym instances.
    class GymWorldGenerator : public Experiment::IInstanceWorld
    {
    public:
        // Constructor does almost nothing. Parameters need to be read in via
        // the readInAllParameters function before the class is useful.
        GymWorldGenerator();
        ~GymWorldGenerator();

        // Inherited from Experiment::IAgentGenerator
        bool loadConfigurationFromYaml(YAML::Node node);
        // After reading in the parameters, instantiate the testing world. This
        // function makes the assumption that readInAllParameters ran successfully.
        bool createTestWorld();

        // Test a population of individuals, return their fitnesses:
        std::vector<double> getPopulationFitnesses(const std::vector<Experiment::IAgent*>& population);

    protected:
        bool m_initialized;
        GYM_TYPE m_gymType;
        YAML::Node m_config;

        // World-specific values we'll read from the config:
        std::string m_environmentName;
        int m_observationSize;
        int m_actionSize;
        bool m_useArgMax;  // Discrete uses arg-max, Continuous doesn't.
        unsigned int m_stepsPerRun;  // How many steps each sim run is allowed before being forceably stopped.
        // Grading values we'll read from config:
        int m_runsPerIndividual;
        int m_runsPerIndividualToGrade;
        FITNESS_COLLAPSE_FUNCTION m_fitnessCollapseFunction;

        // We will only need one of these three:
        CGym::AbstractGymEnv_1D* m_CGym1D;
        CGym::AbstractGymEnv_2D* m_CGym2D;

        // Collapse fitnesses down to a single value. This vector will be sorted as
        // part of the process:
        double applyFitnessCollapseFunction(std::vector<double>& fitnesses);

        // We need to implement our own way of gathering fitnesses from CGym environments:
        std::vector<double> getCGym1DPopulationFitnesses(const std::vector<Experiment::IAgent*>& population);
        double getCGym1DSingleFitness(Experiment::IAgent* agent);
        
        std::vector<double> getCGym2DPopulationFitnesses(const std::vector<Experiment::IAgent*>& population);
        double getCGym2DSingleFitness(Experiment::IAgent* agent);
    private:
        virtual std::string classname() { return "GymWorldGenerator"; }
    };
}

