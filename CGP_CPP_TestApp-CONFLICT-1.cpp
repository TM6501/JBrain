// CGP_CPP_TestApp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <unordered_map>
#include <chrono>
#include <thread>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include "CGPFunctions.h"
#include "StandardCGPIndividual.h"
#include "StandardCGPIndividual_impl.h"
#include "CGPFFANNIndividual.h"
#include "GymTester.h"
#include "GeneralCGPSolver.h"
#include "GymTester.h"
#include "Python.h"
#include "DungeonRoomEnv.h"
#include "SimpleRepeat.h"
#include "IAgent.h"

#include "json.hpp"
using json = nlohmann::json;

#include "yaml-cpp/yaml.h"

#include "Experiment.h"
#include "MultiExperiment.h"
#include "JBrainFactory.h"
#include "JBrain.h"
#include "iostream"
#include "fstream"

#include "GymSageRunner.h"
// 5 - Full logging, all steps
// 0 - No logging
unsigned int DEBUG_LEVEL = 0;

// Still need to run this test of the write/load from JSON:

int testWriteAndLoadOfJSON()
{
    Experiment::Experiment testExperiment("testExperiment.yaml");    
    testExperiment.runExperiment();

    json jOut;

    Experiment::IAgent* finalAgent = testExperiment.getFinalAgent();
    finalAgent->writeSelfToJson(jOut);

    std::cout << "Running test on final agent..." << std::endl;
    double finalScore = testExperiment.testSingleAgent(finalAgent);
    std::cout << "Final score: " << finalScore << std::endl;

    std::cout << "Writing final agent out to file..." << std::endl;

    std::ofstream outFile("testJsonOut.json");
    outFile << jOut << std::endl;
    outFile.close();

    // Once we ask for the final agent from the experiment, we become
    // responsible for its cleanup:
    delete finalAgent;
    finalAgent = nullptr;

    std::cout << "Reading agent back in from file..." << std::endl;
    finalAgent = testExperiment.getAgentFromJson("testJsonOut.json");
        
    double finalScore2 = testExperiment.testSingleAgent(finalAgent);
    std::cout << "Final score: " << finalScore2 << std::endl;    

    return 0;
}

// Linux likes "strcasecmp", Windows complains if it isn't "_stricmp". Since
// this main is so simple, we'll just make a new one for each OS:
#ifdef __linux__

int main(int argc, char** argv)
{
    // Full argument set for multi-run:
    // app.exe <"multiSetup"|"multiRun"> <fileName|folderName> [<numProc>]
    // For single run:
    // app.exe single yamlFilename

    // app.exe multiSetup filename
    if (argc == 3 && strcasecmp(argv[1], "multiSetup") == 0)
    {
        // Separate the multi-experiment-yaml to many single-experiment yaml:
        Experiment::MultiExperiment testExperiment;

        std::string yamlDir = testExperiment.setupExperiments(std::string(argv[2]));

        std::cout << "Multi-experiment setup." << std::endl;
        std::cout << "To run the full folder of experiments, use these arguments next time: " << std::endl;
        std::cout << "<executableName> multiRun \"" << yamlDir << "\" <numParallelProc>" << std::endl;
    }
    // app.exe multiRun folderName <numProc>
    else if (argc == 4 && strcasecmp(argv[1], "multiRun") == 0)
    {
        Experiment::MultiExperiment testExperiment;
        unsigned int numProc = static_cast<unsigned int>(atoi(argv[3]));

        testExperiment.runDirectoryOfExperiments(argv[2], numProc);
    }

    // progName.exe single fileName
    else if (argc == 3 && strcasecmp(argv[1], "single") == 0)
    {
        Experiment::Experiment singleExperiment(argv[2]);
        singleExperiment.runExperiment();
    }
    else
    {
        std::cout << "Usage: \"<appName> <multi|single> <yamlConfig> [numProc if multiExp]\"" << std::endl;
        return -1;
    }
}

#elif _WIN32
/*
int main(int argc, char** argv)
{   
    // Full argument set for multi-run:
    // app.exe <"multiSetup"|"multiRun"> <fileName|folderName> [<numProc>]
    // For single run:
    // app.exe single yamlFilename
        
    // app.exe multiSetup filename
    if (argc == 3 && _stricmp(argv[1], "multiSetup") == 0)
    {
        // Separate the multi-experiment-yaml to many single-experiment yaml:
        Experiment::MultiExperiment testExperiment;

        std::string yamlDir = testExperiment.setupExperiments(std::string(argv[2]));

        std::cout << "Multi-experiment setup." << std::endl;
        std::cout << "To run the full folder of experiments, use these arguments next time: " << std::endl;
        std::cout << "<executableName> multiRun \"" << yamlDir << "\" <numParallelProc>" << std::endl;
    }
    // app.exe multiRun folderName <numProc>
    else if (argc == 4 && _stricmp(argv[1], "multiRun") == 0)
    {
        Experiment::MultiExperiment testExperiment;
        unsigned int numProc = static_cast<unsigned int>(atoi(argv[3]));

        testExperiment.runDirectoryOfExperiments(argv[2], numProc);
    }

    // progName.exe single fileName
    else if (argc == 3 && _stricmp(argv[1], "single") == 0)
    {
        Experiment::Experiment singleExperiment(argv[2]);
        singleExperiment.runExperiment();
    }
    else
    {
        std::cout << "Usage: \"<appName> <multi|single> <yamlConfig> [numProc if multiExp]\"" << std::endl;
        return -1;
    }
}
*/

int randInt(int low, int high)
{
    return (rand() % (high - low + 1)) + low;
}

float randFloat(float low, float high)
{
    return low + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (high - low)));
}

// Test storing a full brain DNA as YAML:
int testStoringBrainAsYAML()
{
    // Standard random variables are used a lot. May eventually decide to use the newer,
    // more sophisticated random number generators.
    srand(static_cast<unsigned int>(time(NULL)));

    YAML::Node root;
    YAML::Node dendrites;
    std::string denInfo = "dendriteInfo";
    root[denInfo]["maxLength"] = 22;
    root["dendriteInfo"]["maxCount"] = 5;
    root["dendriteInfo"]["weights"].push_back(-0.9);
    root["dendriteInfo"]["weights"].push_back(0.9);
    root["dendriteInfo"]["weights"].push_back(0.5);

    YAML::Node neuron1, neuron2;
    neuron1["X"] = 1.2;
    neuron1["Y"] = 0.7;
    neuron1["Z"] = -0.9;
    neuron2["X"] = -1.2;
    neuron2["Y"] = -0.7;
    neuron2["Z"] = 0.9;
    root["neurons"].push_back(neuron1);
    root["neurons"].push_back(neuron2);

    std::ofstream fout("testOut.yaml");
    fout << root;
    fout.close();
}

int testBrainCreationAndPrintOut()
{    
    // JBrain::JBrainFactory* factory = JBrain::JBrainFactory::getInstance();
    // factory->initialize("brainConfig.yaml");
    // JBrain::JBrain* testBrain = factory->getRandomBrain();
 
    // testBrain->writeSelfHumanReadable(std::cout);

    // *********** Test writing brain to JSON *********************
    
    Experiment::GymSageRunner* sageRunner = Experiment::GymSageRunner::getInstance();
    sageRunner->initialize(true, 4, 1);

    std::cout << "If we made it here, the sage runner properly initialized." << std::endl;

    std::vector<double> obs;
    std::vector<double> act;

    if (sageRunner->reset())
    {
        obs = sageRunner->getRecentObs();
        act = sageRunner->getRecentAgentAction();

        std::cout << "Observation: ";
        for (int i = 0; i < obs.size(); ++i)
            std::cout << obs[i] << " ";
        std::cout << std::endl;

        std::cout << "Agent's action: ";
        for (int i = 0; i < act.size(); ++i)
            std::cout << act[i] << " ";
        std::cout << std::endl;
    }
    else
    {
        std::cout << "Sage runner reset failed." << std::endl;
    }

    bool done = false;
    float reward = 0.0;

    while (!done)
    {
        // sageRunner->render();
        
        // Take the action recommended by the sage:
        if (!sageRunner->step(act))
        {
            std::cout << "Step function failed." << std::endl;
            return -1;
        }
        std::cout << "2." << std::endl;

        // Gather information:
        sageRunner->getAllStatus(obs, act, done, reward);

        std::cout << "3." << std::endl;

        std::cout << "Most recent observaction:";
        for (int i = 0; i < obs.size(); ++i)
            std::cout << obs[i] << " ";
        std::cout << std::endl;

        std::cout << "4." << std::endl;

        std::cout << "Recommended action: " << act[0] << std::endl;
        std::cout << "Current reward: " << reward << std::endl;
        std::cout << "Environment done: " << done << std::endl;
    }

    std::cout << "Returning." << std::endl;
    return 0;
}

int main(int argc, char** argv)
{
    return testBrainCreationAndPrintOut();
}

#endif



