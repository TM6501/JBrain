// CGP_CPP_TestApp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <iterator>
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
#include <algorithm>
#include <cmath>

// Needed for experiment data directory manipulation:
#include <iomanip>
#include <filesystem>
#include <ctime>
#include <chrono>
#include <sstream>

#include "CGPFunctions.h"
#include "StandardCGPIndividual.h"
#include "StandardCGPIndividual_impl.h"
#include "CGPFFANNIndividual.h"
#include "GymTester.h"
#include "GeneralCGPSolver.h"
#include "GymTester.h"
#include "Python.h"
// #include "DungeonRoomEnv.h"
// #include "SimpleRepeat.h"
#include "IAgent.h"

#include "JsonLib/json.hpp"
using json = nlohmann::json;

#include "yaml-cpp/yaml.h"

#include "Experiment.h"
#include "MultiExperiment.h"
#include "JBrainFactory.h"
#include "ObservationProcessor.h"
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

// Calculate the slope of a linear trendline through the scores
// to use as an evaluation metric.
// Equation taken from:
// https://math.stackexchange.com/questions/204020/what-is-the-equation-used-to-calculate-a-linear-trendline
double trendlineSlope(const std::vector<double>& scores)
{
    double sumXY = 0.0;
    double sumX = 0.0;
    double sumY = 0.0;
    double sumXsquared = 0.0;

    // For the purposes of this equation, X is the sample number
    // (starting from 1) and Y is the score.
    double x;
    for (unsigned int i = 0; i < scores.size(); ++i)
    {
        x = static_cast<double>(i + 1.0);
        sumXY += (x * scores[i]);
        sumX += x;
        sumY += scores[i];
        sumXsquared += (x * x);
    }
    double count = static_cast<double>(scores.size());
    double slope = ((count * sumXY) - (sumX * sumY)) /
        ((count * sumXsquared) - (sumX * sumX));

    return slope;
}

double testOneBrain(JBrain::JBrain* brain, Experiment::GymSageRunner* sageRunner,
    const std::string& dataDir, const unsigned int& trainingRuns,
    const unsigned int& testingRuns, const double& sageMatchReward,
    const double& trendlineReward, const double& maxMinDiffReward,
    const double& brainCrashReward)
{
    // First, have the brain write itself out to json:
    json jOut;
    brain->writeSelfToJson(jOut);
    std::string jsonName = dataDir + brain->getName() + "_initial.json";
    std::ofstream outFile(jsonName.c_str());
    // outFile << std::setw(2) << jOut << std::endl; // Human readable
    outFile << jOut << std::endl; // Save space
    outFile.close();

    // Initialize our CSV output:
    brain->initializeCSVOutputFile(dataDir);

    // Let it run and hopefully learn:    
    std::vector<double> rewards;
    std::vector<double> testRewards;
    std::vector<double> obs;
    std::vector<double> sageAct;
    int sageChoice;
    std::vector<double> brainAct;
    std::vector<double> brainActIn {0.0};  // The action passed to the step function.
    bool done = false;
    float reward = 0.0;
    static double MAX_REWARD = 500.0; // Hardcoded for cartpole-v1
    bool firstObservation;
    unsigned int testSageMatchChances = 0;
    unsigned int testSageMatches = 0;
    
    for (unsigned int i = 0; i < trainingRuns + testingRuns; ++i)
    {
        done = false;
        sageRunner->reset();
        sageRunner->getAllStatus(obs, sageAct, done, reward);
        
        firstObservation = true;
        while (!done)
        {
            // Simple environment, sage choice is 0 or 1:
            sageChoice = 0;
            if (sageAct[0] > 0.5)
                sageChoice = 1;

            // Let the brain decide what to do:
            brainAct = brain->processInput(obs, sageChoice, firstObservation);
            firstObservation = false;
         
            // The environment expects a single element double vector with
            // 0 and 1 being the only acceptable values:
            brainActIn[0] = 0; // Assume 0, change if needed:
            if (brainAct[1] > brainAct[0])
                brainActIn[0] = 1.0;

            // Check sage-matching:
            if (i >= trainingRuns)
            {
                ++testSageMatchChances;
                // Both 1:
                if (sageChoice == 1 && brainAct[1] > brainAct[0])
                    ++testSageMatches;

                // Both 0:
                else if (sageChoice == 0 && brainAct[1] <= brainAct[0])
                    ++testSageMatches;
            }

            sageRunner->step(brainActIn);

            // Gather information:
            sageRunner->getAllStatus(obs, sageAct, done, reward);
        }

        // Tell the brain to write out how well it did:
        brain->writeLineToCSVOutputFile(reward);

        if (i >= trainingRuns)
            testRewards.push_back(reward);

        rewards.push_back(reward);
    }

    // Have the brain close out its csv:
    brain->closeCSVOutputFile();

    // Write out what the brain looked like in the end:
    json jOut2;
    brain->writeSelfToJson(jOut2);
    std::string jsonName2 = dataDir + brain->getName() + "_final.json";
    std::ofstream outFile2(jsonName2.c_str());
    // outFile << std::setw(2) << jOut << std::endl; // Human readable
    outFile2 << jOut2 << std::endl; // Save space
    outFile2.close();

    // Calculate each reward:
    double fullReward = trendlineSlope(rewards) * trendlineReward;

    // Get the min and max of the test rewards:
    double minTest = testRewards[0];
    double maxTest = testRewards[0];
    for (double& rew : testRewards)
    {
        minTest = fmin(minTest, rew);
        maxTest = fmax(maxTest, rew);
    }
    double testDiff = (maxTest / MAX_REWARD) - (minTest / MAX_REWARD);
    fullReward += testDiff * maxMinDiffReward;

    double sageMatchPercent = static_cast<double>(testSageMatches) /
        static_cast<double>(testSageMatchChances);

    fullReward += sageMatchPercent * sageMatchReward;

    // Check for brain crash:
    if (brain->getNeuronCount() == 0)
        fullReward += brainCrashReward;

    return fullReward;
    /*
    // Final score is the average score during test time:
    double finalReward;
    double averageReward = std::accumulate(rewards.begin(), rewards.end(), 0.0) / rewards.size();

    // Not a true variance/stdDev calculation, but need to punish having wild
    // variety. The goal is stability:
    finalReward = averageReward;
    for (unsigned int i = 0; i < rewards.size(); ++i)
        finalReward -= fabs(averageReward - rewards[0]) / testingRuns;

    // Lose points if the brain died:
    if (brain->getNeuronCount() == 0)
        finalReward -= 500.0;

    return finalReward;
    */
}

void getInfoFromInProgressExperiment(const std::string& directory,
    int& previousEpoch, JBrain::JBrain*& parent)
{
    previousEpoch = 0;
    parent = nullptr;

    // Open the CSV from the data directory for reading:
    std::ifstream csv(directory + "\\experimentSummary.csv");
    if (!csv)
    {
        std::cout << "experimentSummary.csv not found in " << directory << std::endl;
        return;
    }

    std::string line;
    std::string prevLine;
    std::string cell;

    // Read in lines until we reach the bottom of the csv:
    while (std::getline(csv, line)) { prevLine = line; }

    // The final line should now be in "line":
    line = prevLine;

    // Get what we need from the line:
    std::istringstream lineStream(line);
    std::string brainFileName = "";
    unsigned int idx = 0;

    // This function will need to be changed if the column index of 
    // the epoch (0) or the best brain (2) changes:
    while(std::getline(lineStream, cell, ','))
    {
        if (idx == 0)
            previousEpoch = std::stoi(cell);

        else if (idx == 2)
        {
            std::ifstream jsonInFile(brainFileName = directory + "\\" + cell + "_final.json");
            auto brainJson = json::parse(jsonInFile);
            jsonInFile.close();

            parent = JBrain::JBrain::getBrainFromJson(brainJson);
            break;
        }
        ++idx;
    }
}

int testFullExperiment(std::string yamlFileName)
{
    YAML::Node fullConfig = YAML::LoadFile(yamlFileName);
    YAML::Node expConfig;
    if (fullConfig["Experiment"])
        expConfig = fullConfig["Experiment"];
    else
    {
        std::cout << "Experiment config not defined." << std::endl;
        return -1;
    }

    std::string dataDir;

    // If an experiment-in-progress was provided, use it:
    bool expInProgress = false;
    JBrain::JBrain* parent = nullptr;
    int prevEpoch = -1;
    std::string expDir = expConfig["ExperimentDirectory"].as<std::string>();
    if (expDir != "None")
    {
        expInProgress = true;
        dataDir = expDir;

        // Load up the new yaml and change the yamlFileName:
        yamlFileName = dataDir + "experimentConfigFile.yaml";
        fullConfig = YAML::LoadFile(yamlFileName);
        expConfig = fullConfig["Experiment"];

        getInfoFromInProgressExperiment(dataDir, prevEpoch, parent);
    }
    else // Need to build the directory:
    {
        std::string mainDir = expConfig["MainDataDirectory"].as<std::string>();
        std::string expName = expConfig["ExperimentName"].as<std::string>();

        // I can't believe C++ has its head up its ass so far that this is required
        // to get the friggen local time:
        auto now = std::chrono::system_clock::now();
        time_t t = std::chrono::system_clock::to_time_t(now);
        std::tm tm{};
#if defined(__unix__)
        localtime_r(&t, &tm);
#else
        localtime_s(&tm, &t);
#endif

        // Get the date-time as a string for a unique experiment directory name:
        char buffer[80];
        strftime(buffer, sizeof(buffer), "_%Y-%m-%d_%H-%M-%S", &tm);

        // Combine them into a single directory name:
        dataDir = mainDir + expName + std::string(buffer) + "\\";
        std::filesystem::create_directories(dataDir);

        // Copy our configuration file into the experiment directory:
        std::filesystem::copy_file(yamlFileName, dataDir + "experimentConfigFile.yaml");
    }

    // Create a full experiment CSV, even if we have nothing to write right now:
    std::string csvFName = dataDir + "experimentSummary.csv";
    std::ofstream csvOut(csvFName.c_str(), std::ios_base::app);  // Append if it exists already.

    // Write the header if we're just getting started:
    if (!expInProgress)
        csvOut << "epoch,bestScore,bestBrain,bestParent,bestNeuronCount,trainingSeconds" << std::endl;

    csvOut.close();  // Reopen for each write since they are infrequent.

    Experiment::GymSageRunner* sageRunner = Experiment::GymSageRunner::getInstance();
    // useArgMax MUST be false:
    sageRunner->initialize(false, 4, 1);
       
    std::cout << "Initializing the brain factory..." << std::endl;
    JBrain::JBrainFactory* factory = JBrain::JBrainFactory::getInstance();    
    factory->initialize(yamlFileName);
    
    unsigned int startPop = expConfig["StartingPopulationSize"].as<unsigned int>();
    std::vector<JBrain::JBrain*> population;

    if (expInProgress)
    {
        std::cout << "Generating population from previous best brain..." << std::endl;
        // New brain number to start from:
        std::string brainName = parent->getName();
        unsigned int nextBrainNumber = static_cast<unsigned int>(
            std::stoi(brainName.substr(1)));

        factory->setNextBrainNumber(nextBrainNumber + 1);
        population = factory->getFullMutatedPopulation(parent);
    }
    else
    {
        std::cout << "Getting " << startPop << " random sample brains..." << std::endl;
        for (unsigned int i = 0; i < startPop; ++i)
            population.push_back(factory->getRandomBrain());
    }

    unsigned int trainingEpochs = expConfig["MaximumEpochs"].as<unsigned int>();
    std::vector<double> allRewards;
    double endReward = expConfig["MaximumReward"].as<double>();
    double minReward = expConfig["PopulationInfusionReward"].as<double>();
    unsigned int populationAddSize = expConfig["PopulationInfusionSize"].as<unsigned int>();
    unsigned int trainingRuns = expConfig["TrainTrials"].as<unsigned int>();
    unsigned int testingRuns = expConfig["TestTrials"].as<unsigned int>();

    // Available reward calculations:
    double sageMatchReward = expConfig["Reward_TestPercentSageMatch"].as<double>();
    double trendSlopeReward = expConfig["Reward_AllTrendlineSlope"].as<double>();
    double variationReward = expConfig["Reward_TestMaxMinPercentDiff"].as<double>();
    double brainCrashReward = expConfig["Reward_PenaltyForBrainCrash"].as<double>();

    double maxReward;
    int maxIndex;
    std::vector<double> obs;
    std::vector<double> sageAct;
    std::vector<double> brainAct;
    bool done = false;
    float reward = 0.0;
    for (unsigned int i = static_cast<unsigned int>(prevEpoch + 1); i < trainingEpochs; ++i)
    {
        auto start = std::chrono::high_resolution_clock::now();
        
        std::cout << "Beginning training on epoch " << i + 1 << " / " << trainingEpochs;
        allRewards.clear();

        for (auto brain : population)
        {
            allRewards.push_back(testOneBrain(brain, sageRunner, dataDir,
                trainingRuns, testingRuns, sageMatchReward, trendSlopeReward,
                variationReward, brainCrashReward));
            std::cout << ".";
        }

        maxReward = -1.0;
        maxIndex = -1;
        for (unsigned int j = 0; j < allRewards.size(); ++j)
        {
            if (allRewards[j] > maxReward)
            {
                maxReward = allRewards[j];
                maxIndex = j;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
        
        // Gather the best:
        parent = population[maxIndex];

        std::cout << " Max reward: " << maxReward << " from brain " <<
            parent->getName() << ". Training took " << seconds 
            << " seconds." << std::endl;
        
        // Write line to our output CSV file
        // epoch, best score, best brain, best parent, best brain's neuron count,
        // training seconds
        csvOut.open(csvFName.c_str(), std::ios_base::app);
        csvOut << i << ","
            << maxReward << ","
            << parent->getName() << ","
            << parent->getParentName() << ","
            << parent->getNeuronCount() << ","
            << seconds << std::endl;
        csvOut.close();

        // Delete all but the best:
        population.erase(population.begin() + maxIndex);
        for (auto brain : population)
            delete brain;

        // Create the new generation:
        population.clear();
        population = factory->getFullMutatedPopulation(parent);

        // If they're doing poorly, add in some fresh blood:
        if (maxReward < minReward)
        {
            std::cout << "Poor overall performance. " << maxReward
                << " < " << minReward << ". Adding " << populationAddSize
                << " random individuals to the population." << std::endl;

            for (unsigned int i = 0; i < populationAddSize; ++i)
                population.push_back(factory->getRandomBrain());
        }

        // End the experiment early:
        if (maxReward > endReward)
        {
            std::cout << "Reward greater than " << endReward << " achieved. Ending." << std::endl;
            i = trainingEpochs + 5;
        }
    }

    // Clean up:
    for (auto brain : population)
        delete brain;

    return 0;
}

int testGymSageRunner_lunarLander()
{
  Experiment::GymSageRunner* sageRunner = Experiment::GymSageRunner::getInstance();
  sageRunner->initialize(false, 8, 1, "LunarLander-v3", "ppo_sage_lunarLanderDiscrete", "human");

  std::cout << "Testing string arguments." << std::endl;

  sageRunner->testStringArguments("StringArg1", "StringArg4_JustKidding_2");

  std::cout << "Back from testStringArguments." << std::endl;
  
  // std::cout << "If we made it here, the sage runner properly initialized." << std::endl;

  std::vector<double> obs;
  std::vector<double> act;

  bool done = false;
  float reward = 0.0;
  std::vector<float> allRewards;
  int timesToRun = 5;

  std::cout << "Beginning " << timesToRun << " runs of LunarLander-v3." << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < timesToRun; ++i)
  {
    sageRunner->reset();
    sageRunner->getAllStatus(obs, act, done, reward);

    while (!done)
    {
      // sageRunner->render();

      // Take the action recommended by the sage:
      if (!sageRunner->step(act))
      {
        std::cout << "Step function failed." << std::endl;
        return -1;
      }

      // Gather information:
      sageRunner->getAllStatus(obs, act, done, reward);

      /*
      std::cout << "Next action::";
      for (int i = 0; i < act.size(); ++i)
          std::cout << act[i] << ", ";
      std::cout << std::endl;
      */

      // std::cout << "TestApp recommended action: " << act[0] << std::endl;
      // std::cout << "Current reward: " << reward << std::endl;
      // std::cout << "Environment done: " << done << std::endl;
    }
    allRewards.push_back(reward);
  }
  auto end = std::chrono::high_resolution_clock::now();

  float averageReward = static_cast<float>(std::accumulate(allRewards.begin(),
    allRewards.end(), 0) / timesToRun);

  std::cout << timesToRun << " runs through LunarLander-v3 took " <<
    std::chrono::duration_cast<std::chrono::seconds>(end - start).count() <<
    " seconds. The average run got a reward of " << averageReward << "." <<
    std::endl;

  return 0;
}

int testGymSageRunner_cartPole()
{
    Experiment::GymSageRunner* sageRunner = Experiment::GymSageRunner::getInstance();
    sageRunner->initialize(false, 4, 1, "CartPole-v1", "ppo_sage", "human");
    sageRunner->reset();
    return 0;
    // std::cout << "If we made it here, the sage runner properly initialized." << std::endl;

    std::vector<double> obs;
    std::vector<double> act;
   
    bool done = false;
    float reward = 0.0;
    std::vector<float> allRewards;
    int timesToRun = 5;

    std::cout << "Beginning " << timesToRun << " runs of CartPole-v1." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < timesToRun; ++i)
    {
        sageRunner->reset();
        sageRunner->getAllStatus(obs, act, done, reward);
         

        while (!done)
        {
            sageRunner->render();

            // Take the action recommended by the sage:
            if (!sageRunner->step(act))
            {
                std::cout << "Step function failed." << std::endl;
                return -1;
            }

            // Gather information:
            sageRunner->getAllStatus(obs, act, done, reward);

            /*
            std::cout << "Next action::";
            for (int i = 0; i < act.size(); ++i)
                std::cout << act[i] << ", ";
            std::cout << std::endl;
            */

            // std::cout << "TestApp recommended action: " << act[0] << std::endl;
            // std::cout << "Current reward: " << reward << std::endl;
            // std::cout << "Environment done: " << done << std::endl;
        }
        allRewards.push_back(reward);
    }
    auto end = std::chrono::high_resolution_clock::now();

    float averageReward = static_cast<float>(std::accumulate(allRewards.begin(),
        allRewards.end(), 0) / timesToRun);

    std::cout << timesToRun << " runs through CartPole-v1 took " <<
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count() <<
        " seconds. The average run took " << averageReward << " steps to complete." <<
        std::endl;

    return 0;
}

int testMutation()
{
    std::cout << "Initializing the brain factory..." << std::endl;
    JBrain::JBrainFactory* factory = JBrain::JBrainFactory::getInstance();
    factory->initialize("brainConfig.yaml");

    std::cout << "Getting a random sample brain..." << std::endl;
    JBrain::JBrain* testBrain = factory->getRandomBrain();

    json jOut;
    testBrain->writeSelfToJson(jOut);
    std::ofstream outFile("originalBrain.json");
    outFile << std::setw(2) << jOut << std::endl;
    outFile.close();
    
    std::vector<JBrain::JBrain*> nextGen = factory->getFullMutatedPopulation(testBrain);

    std::cout << nextGen.size() << " additional brains created." << std::endl;

    for (unsigned int i = 0; i < nextGen.size(); ++i)
    {
        json jOut;
        nextGen[i]->writeSelfToJson(jOut);
        std::stringstream filename;
        filename << i << ".json";
        std::ofstream outFile(filename.str());
        outFile << std::setw(2) << jOut << std::endl;
        outFile.close();
        delete nextGen[i];
        nextGen[i] = nullptr;
    }

    return 0;
}

#include <random>
#include <numeric>
int testProbabilitySelection()
{
  std::mt19937_64 gen(std::random_device{}());
  std::vector<unsigned int> chances{ 1, 1, 1, 1, 1 };
  std::vector<unsigned int> selections(5);
  std::iota(selections.begin(), selections.end(), 0);

  std::discrete_distribution<std::size_t> dist{ chances.begin(), chances.end() };

  std::cout << "Selecting 1000 values: ";
  for (unsigned int i = 0; i < 100; ++i)
    std::cout << selections[dist(gen)] << ", ";
  std::cout << selections[dist(gen)] << std::endl;

  return 0;
}

#include "ObservationProcessor.h"
#include <cstdlib>
int testObservationProcessor()
{
  std::vector<std::vector<double> > obsRanges{
    std::vector<double>{-2.4f, 2.4f},
    std::vector<double>{-1.0f, 1.0f},
    std::vector<double>{-0.3f, 0.3f},
    std::vector<double>{-1.0f, 1.0f}
  };

  std::vector<std::vector<double> > actRanges{
    std::vector<double>{-1.0, 1.0},
    std::vector<double>{-10.0, 10.0}
  };

  JBrain::ObservationProcessor obProc_buckets(
    CGP::INPUT_PREPROCESSING::BUCKETS,
    4,
    obsRanges,
    { 10,10,10,10 },
    2,
    actRanges,
    { 100, 100 });    

  std::vector<double> observation;
  std::vector<unsigned int> action;
  double x;
  unsigned int singleAct;
  std::vector<double> stdOut;
  std::vector<double> bucketOut;
  std::vector<std::vector<double> > bucketSep;
  std::vector<unsigned int> bucketSepSimplified;
  std::vector<double> actOut;
  std::vector<unsigned int> outBucketsSimplified;
  while (true)
  {
    observation.clear();
    std::cout << "4 input values: ";
    for (unsigned int i = 0; i < 4; ++i)
    {
      std::cin >> x;
      observation.push_back(x);
    }

    bucketOut = obProc_buckets.processInput(observation);

    std::cout << "Input vector: ";
    for (double elem : observation)
      std::cout << elem << " ";
    std::cout << std::endl;
        
    std::cout << "Bucket processing: ";
    for (double elem : bucketOut)
      std::cout << elem << " ";
    std::cout << std::endl;

    std::cout << "Bucket separated: " << std::endl;
    bucketSep = obProc_buckets.getSeparatedInputs();
    for (auto topLevel : bucketSep)
    {
      std::cout << "\t";
      for (auto singleVal : topLevel)
        std::cout << singleVal << " ";
      std::cout << std::endl;
    }
    
    std::cout << "Bucket separated simplified: ";
    bucketSepSimplified = obProc_buckets.getSeparatedInputs_simplified();
    for (auto val : bucketSepSimplified)
      std::cout << val << " ";
    std::cout << std::endl;


    std::cout << "2 output value buckets: ";
    action.clear();
    for (unsigned int i = 0; i < 2; ++i)
    {
      std::cin >> singleAct;
      action.push_back(singleAct);
    }

    actOut = obProc_buckets.processOutput(action);
    std::cout << "Action output: ";
    for (auto act : actOut)
      std::cout << act << "  ";
    std::cout << std::endl;

    outBucketsSimplified = obProc_buckets.getSeparatedOutputs_simplified(actOut);
    std::cout << "Action buckets: ";
    for (auto actBuck : outBucketsSimplified)
      std::cout << actBuck << " ";
    std::cout << std::endl;
  }

  return 0;
}

bool getEventHappened(double probability)
{
  // Random device and distribution don't need to be
  // recreated every time:
  static std::random_device rd;
  static std::mt19937_64 eng(rd());
  static std::uniform_real_distribution<double> distr(0.0, 1.0);

  if (distr(eng) <= probability)
    return true;
  else
    return false;
}

int testProbability()
{
  unsigned int timesToCheck = 10000000;
  std::vector<double> probsToCheck{ 0.0, 1.0, 0.5, 0.99, 0.01, 0.0001, 0.999, 0.25, 0.75 };

  unsigned int happenedCount;
  for (auto& prob : probsToCheck)
  {
    happenedCount = 0;
    for (unsigned int i = 0; i < timesToCheck; ++i)
      if (getEventHappened(prob))
        ++happenedCount;

    std::cout << "Event with probability of " << prob * 100.0
      << "% happened approximately "
      << (static_cast<double>(happenedCount) / static_cast<double>(timesToCheck)) * 100.0
      << "% of the time." << std::endl;
  }

  return 0;
}

double testSnapJBrain_HDCCSV(JBrain::JBrain_Snap* brain, Experiment::GymSageRunner* sageRunner,
  const std::string& dataDir, const unsigned int& testingRuns,
  std::vector<std::vector<double> >* obsTrainData, std::vector<std::vector<double> >* actTrainData)
{
  // Let it run and hopefully learn:    
  std::vector<double> rewards;
  std::vector<double> testRewards;
  std::vector<double> obs;
  std::vector<double> sageAct;
  std::vector<double> brainAct;
  std::vector<double> brainActIn{ 0.0 };  // The action passed to the step function.
  bool done = false;
  float reward = 0.0;
  // static double MAX_REWARD = 500.0; // Hardcoded for cartpole-v1
  static double MAX_REWARD = 200.0; // Hardcoded for LunarLander-v3
  unsigned int testSageMatchChances = 0;
  unsigned int testSageMatches_temp = 0;
  unsigned int testSageMatches_full = 0;
  unsigned int testSageDiffs_temp = 0;
  unsigned int testSageDiffs_full = 0;
  unsigned int noOutputEvents_temp = 0;
  unsigned int noOutputEvents_full = 0;
  unsigned int procCreated_temp = 0;
  unsigned int procDestroyed_temp = 0;
  unsigned int procCreatedOrDestroyed_full = 0;
  unsigned int inputCreated_temp = 0;
  unsigned int inputDestroyed_temp = 0;
  unsigned int inputCreatedOrDestroyed_full = 0;

  // Initialize the brain's output:
  brain->initializeCSVOutputFile(dataDir);

  std::cout << "\t" << brain->getName() << " beginning training...";

  // Let the brain see each training data sample exactly once:
  auto start = std::chrono::high_resolution_clock::now();
  for (unsigned int idx = 0; idx < obsTrainData->size(); ++idx)
  {
    brain->processInput((*obsTrainData)[idx], (*actTrainData)[idx]);
  }
  auto mid = std::chrono::high_resolution_clock::now();
  auto seconds = std::chrono::duration_cast<std::chrono::seconds>(mid - start).count();
  std::cout << " training took " << seconds << " seconds. Beginning testing...";

  // All of the training data will be recorded as a single
  // trial run in the CSV:
  brain->processEndOfTrial(0.0, 0.0, MAX_REWARD);  // 0 reward.

  // Turn off training:
  // std::string constructor explicitly required for some reason.
  // "NONE" was being implicitly converted to a boolean and the
  // wrong function was being called.
  brain->setValueByName("HDCLearnMode", std::string("NONE"));

  for (unsigned int i = 0; i < testingRuns; ++i)
  {
    done = false;
    sageRunner->reset();
    sageRunner->getAllStatus(obs, sageAct, done, reward);

    while (!done)
    {
      // Round to make sure little double-errors don't cause issues:
      // sageChoice = std::round(sageAct[0]);

      // brainAct will have length equal to the action size:
      brainAct = brain->processInput(obs, sageAct);

      // The environment expects a single element double vector with
      // the equivalent int action in the 0 index:
      brainActIn[0] = std::max_element(brainAct.begin(), brainAct.end()) - brainAct.begin();

      // Record how often we had the chance to choose correctly during test runs:
      ++testSageMatchChances;

      sageRunner->step(brainActIn);

      // Gather information:
      sageRunner->getAllStatus(obs, sageAct, done, reward);
    }

    testRewards.push_back(reward);
    
    brain->getFullTrialStatistics(noOutputEvents_temp, testSageMatches_temp, testSageDiffs_temp,
        procCreated_temp, procDestroyed_temp, inputCreated_temp, inputDestroyed_temp);
    noOutputEvents_full += noOutputEvents_temp;
    testSageMatches_full += testSageMatches_temp;
    testSageDiffs_full += testSageDiffs_temp;
    procCreatedOrDestroyed_full += procCreated_temp + procDestroyed_temp;
    inputCreatedOrDestroyed_full += inputCreated_temp + inputDestroyed_temp;    

    brain->processEndOfTrial(reward, 0.0, MAX_REWARD);
  
    rewards.push_back(reward);
  }

  brain->closeCSVOutputFile();

  json jOut;
  brain->writeSelfToJson(jOut);
  std::string jsonName = dataDir + brain->getName() + "_final.json";
  std::ofstream outFile(jsonName.c_str());
  // outFile << std::setw(2) << jOut << std::endl; // Human readable
  outFile << jOut << std::endl; // Save space
  outFile.close();

  // Full reward calculation:  
  auto minMax = std::minmax_element(testRewards.begin(), testRewards.end());
  double minScore = *minMax.first;
  double maxScore = *minMax.second;
  double avgScorePercent = (std::accumulate(testRewards.begin(), testRewards.end(), 0.0) / static_cast<double>(testRewards.size())) / MAX_REWARD;
  double scoreDiff = maxScore - minScore;
  double scoreDiffPercent = scoreDiff / MAX_REWARD;
  double percGoodOutput = static_cast<double>(testSageMatches_full) / static_cast<double>(testSageMatchChances);
  double percBadOutput = static_cast<double>(testSageDiffs_full) / static_cast<double>(testSageMatchChances);
  double percBestScore = maxScore / MAX_REWARD;  // Assumes for now that min-possible is zero. 
  double percWorstScore = scoreDiff / MAX_REWARD;  // Change later if more general % is needed.
  double procChurn = static_cast<double>(procCreatedOrDestroyed_full) / (2.0 * static_cast<double>(testSageMatchChances));
  double inputChurn = static_cast<double>(inputCreatedOrDestroyed_full) / (2.0 * static_cast<double>(testSageMatchChances));
  // double rewardTrendlineSlope = trendlineSlope(rewards);

  // For HDC mode, we just need how well we did during testing:
  double fullReward = percBestScore + avgScorePercent;

  auto end = std::chrono::high_resolution_clock::now();
  seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
  std::cout << "testing took " << seconds << " seconds." << std::endl;

  return fullReward;
}

double testSnapJBrain_HDC(JBrain::JBrain_Snap* brain, Experiment::GymSageRunner* sageRunner,
  const std::string& dataDir, const unsigned int& trainingRuns, const unsigned int& testingRuns)
{
  // Let it run and hopefully learn:    
  std::vector<double> rewards;
  std::vector<double> testRewards;
  std::vector<double> obs;
  std::vector<double> sageAct;
  std::vector<double> brainAct;
  std::vector<double> brainActIn{ 0.0 };  // The action passed to the step function.
  bool done = false;
  float reward = 0.0;
  // static double MAX_REWARD = 500.0; // Hardcoded for cartpole-v1
  static double MAX_REWARD = 200.0; // Hardcoded for LunarLander-v3
  unsigned int testSageMatchChances = 0;
  unsigned int testSageMatches_temp = 0;
  unsigned int testSageMatches_full = 0;
  unsigned int testSageDiffs_temp = 0;
  unsigned int testSageDiffs_full = 0;
  unsigned int noOutputEvents_temp = 0;
  unsigned int noOutputEvents_full = 0;
  unsigned int procCreated_temp = 0;
  unsigned int procDestroyed_temp = 0;
  unsigned int procCreatedOrDestroyed_full = 0;
  unsigned int inputCreated_temp = 0;
  unsigned int inputDestroyed_temp = 0;
  unsigned int inputCreatedOrDestroyed_full = 0;

  // Initialize the brain's output:
  brain->initializeCSVOutputFile(dataDir);

  for (unsigned int i = 0; i < trainingRuns + testingRuns; ++i)
  {
    done = false;
    sageRunner->reset();
    sageRunner->getAllStatus(obs, sageAct, done, reward);

    // Turn off training:
    if (i >= trainingRuns)
    {
      // std::string constructor explicitly required for some reason.
      // "NONE" was being implicitly converted to a boolean and the
      // wrong function was being called.
      brain->setValueByName("HDCLearnMode", std::string("NONE"));
    }

    // std::cout << "Trial " << i + 1 << ": ";
    while (!done)
    {
      // Round to make sure little double-errors don't cause issues:
      // sageChoice = std::round(sageAct[0]);      
      
      // brainAct will have length equal to the action size:
      brainAct = brain->processInput(obs, sageAct);
            
      // The environment expects a single element double vector with
      // the equivalent int action in the 0 index:
      brainActIn[0] = std::max_element(brainAct.begin(), brainAct.end()) - brainAct.begin();      

      // Record how often we had the chance to choose correctly during test runs:
      if (i >= trainingRuns)
        ++testSageMatchChances;

      sageRunner->step(brainActIn);

      // Gather information:
      sageRunner->getAllStatus(obs, sageAct, done, reward);
    }

    if (i >= trainingRuns)
    {
      testRewards.push_back(reward);

      brain->getFullTrialStatistics(noOutputEvents_temp, testSageMatches_temp, testSageDiffs_temp,
        procCreated_temp, procDestroyed_temp, inputCreated_temp, inputDestroyed_temp);
      noOutputEvents_full += noOutputEvents_temp;
      testSageMatches_full += testSageMatches_temp;
      testSageDiffs_full += testSageDiffs_temp;
      procCreatedOrDestroyed_full += procCreated_temp + procDestroyed_temp;
      inputCreatedOrDestroyed_full += inputCreated_temp + inputDestroyed_temp;
    }

    brain->processEndOfTrial(reward, 0.0, MAX_REWARD);

    // std::cout << reward << " " << std::endl;   
    rewards.push_back(reward);
  }

  brain->closeCSVOutputFile();

  json jOut;
  brain->writeSelfToJson(jOut);
  std::string jsonName = dataDir + brain->getName() + "_final.json";
  std::ofstream outFile(jsonName.c_str());
  // outFile << std::setw(2) << jOut << std::endl; // Human readable
  outFile << jOut << std::endl; // Save space
  outFile.close();

  // Full reward calculation:  
  auto minMax = std::minmax_element(testRewards.begin(), testRewards.end());
  double minScore = *minMax.first;
  double maxScore = *minMax.second;
  double avgScorePercent = (std::accumulate(testRewards.begin(), testRewards.end(), 0.0) / static_cast<double>(testRewards.size())) / MAX_REWARD;
  double scoreDiff = maxScore - minScore;
  double scoreDiffPercent = scoreDiff / MAX_REWARD;
  double percGoodOutput = static_cast<double>(testSageMatches_full) / static_cast<double>(testSageMatchChances);
  double percBadOutput = static_cast<double>(testSageDiffs_full) / static_cast<double>(testSageMatchChances);
  double percBestScore = maxScore / MAX_REWARD;  // Assumes for now that min-possible is zero. 
  double percWorstScore = scoreDiff / MAX_REWARD;  // Change later if more general % is needed.
  double procChurn = static_cast<double>(procCreatedOrDestroyed_full) / (2.0 * static_cast<double>(testSageMatchChances));
  double inputChurn = static_cast<double>(inputCreatedOrDestroyed_full) / (2.0 * static_cast<double>(testSageMatchChances));
  double rewardTrendlineSlope = trendlineSlope(rewards);

  // For HDC mode, we just need how well we did during testing:
  double fullReward = percBestScore + avgScorePercent - scoreDiffPercent;

  return fullReward;
}

double testSnapJBrain(JBrain::JBrain_Snap* brain, Experiment::GymSageRunner* sageRunner,
  const std::string& dataDir, const unsigned int& trainingRuns, const unsigned int& testingRuns,
  const double& reward_testNoOutput, const double& reward_testPercentGoodOutput,
  const double& reward_testPercentBadOutput, const double& reward_testMaxScorePercent,
  const double& reward_testBestMinusMin, const double& reward_testProcChurn,
  const double& reward_testInputChurn, const double& reward_slopeLinearTrendline)
{
  // Let it run and hopefully learn:    
  std::vector<double> rewards;
  std::vector<double> testRewards;
  std::vector<double> obs;
  std::vector<double> sageAct;
  std::vector<double> brainAct;
  std::vector<double> brainActIn{ 0.0 };  // The action passed to the step function.
  bool done = false;
  float reward = 0.0;
  static double MAX_REWARD = 500.0; // Hardcoded for cartpole-v1
  unsigned int testSageMatchChances = 0;
  unsigned int testSageMatches_temp = 0;
  unsigned int testSageMatches_full = 0;
  unsigned int testSageDiffs_temp = 0;
  unsigned int testSageDiffs_full = 0;
  unsigned int noOutputEvents_temp = 0;
  unsigned int noOutputEvents_full = 0;
  unsigned int procCreated_temp = 0;
  unsigned int procDestroyed_temp = 0;
  unsigned int procCreatedOrDestroyed_full = 0;
  unsigned int inputCreated_temp = 0;
  unsigned int inputDestroyed_temp = 0;
  unsigned int inputCreatedOrDestroyed_full = 0;
  
  // Initialize the brain's output:
  brain->initializeCSVOutputFile(dataDir);
  
  for (unsigned int i = 0; i < trainingRuns + testingRuns; ++i)
  {
    done = false;
    sageRunner->reset();
    sageRunner->getAllStatus(obs, sageAct, done, reward);

    // std::cout << "Trial " << i + 1 << ": ";
    while (!done)
    {
      // Simple environment, sage choice is 0 or 1:
      // sageChoice = 0;
      // if (sageAct[0] > 0.5)
      //  sageChoice = 1;

      // Let the brain decide what to do:
      // std::cout << ".";
      brainAct = brain->processInput(obs, sageAct);      

      // std::cout << "[" << brainAct[0] << ", " << brainAct[1] << "] ";
      // The environment expects a single element double vector with
      // 0 and 1 being the only acceptable values:
      brainActIn[0] = 0; // Assume 0, change if needed:
      if (brainAct[1] > brainAct[0])
        brainActIn[0] = 1.0;

      // Record how often we had the chance to choose correctly during test runs:
      if (i >= trainingRuns)
        ++testSageMatchChances;

      sageRunner->step(brainActIn);

      // Gather information:
      sageRunner->getAllStatus(obs, sageAct, done, reward);
    }
    
    if (i >= trainingRuns)
    {
      testRewards.push_back(reward);

      brain->getFullTrialStatistics(noOutputEvents_temp, testSageMatches_temp, testSageDiffs_temp,
        procCreated_temp, procDestroyed_temp, inputCreated_temp, inputDestroyed_temp);
      noOutputEvents_full += noOutputEvents_temp;
      testSageMatches_full += testSageMatches_temp;
      testSageDiffs_full += testSageDiffs_temp;
      procCreatedOrDestroyed_full += procCreated_temp + procDestroyed_temp;
      inputCreatedOrDestroyed_full += inputCreated_temp + inputDestroyed_temp;
    }

    brain->processEndOfTrial(reward, 0.0, MAX_REWARD);
    
    // std::cout << reward << " " << std::endl;   
    rewards.push_back(reward);
  }

  brain->closeCSVOutputFile();

  json jOut;
  brain->writeSelfToJson(jOut);
  std::string jsonName = dataDir + brain->getName() + "_final.json";
  std::ofstream outFile(jsonName.c_str());
  // outFile << std::setw(2) << jOut << std::endl; // Human readable
  outFile << jOut << std::endl; // Save space
  outFile.close();

  // Full reward calculation:  
  auto minMax = std::minmax_element(testRewards.begin(), testRewards.end());
  double minScore = *minMax.first;
  double maxScore = *minMax.second;
  double scoreDiff = maxScore - minScore;
  double percGoodOutput = static_cast<double>(testSageMatches_full) / static_cast<double>(testSageMatchChances);
  double percBadOutput = static_cast<double>(testSageDiffs_full) / static_cast<double>(testSageMatchChances);
  double percBestScore = maxScore / MAX_REWARD;  // Assumes for now that min-possible is zero. 
  double percWorstScore = scoreDiff / MAX_REWARD;  // Change later if more general % is needed.
  double procChurn = static_cast<double>(procCreatedOrDestroyed_full) / (2.0 * static_cast<double>(testSageMatchChances));
  double inputChurn = static_cast<double>(inputCreatedOrDestroyed_full) / (2.0 * static_cast<double>(testSageMatchChances));
  double rewardTrendlineSlope = trendlineSlope(rewards);

  // fullReward = %GoodOutput + %BadOutput + (bestScore / MaxScore) + (worstScore / MaxScore) +
  //              noOutputEventCount + procNeuronChurn + inputNeuronChurn.
  // With each value having a multiplier, most likely negative for the bad characteristics:
  double fullReward =
    (percGoodOutput * reward_testPercentGoodOutput) +
    (percBadOutput * reward_testPercentBadOutput) +
    (percBestScore * reward_testMaxScorePercent) +
    (percWorstScore * reward_testBestMinusMin) +
    (static_cast<double>(noOutputEvents_full) * reward_testNoOutput) +
    (procChurn * reward_testProcChurn) +
    (inputChurn * reward_testInputChurn) +
    (rewardTrendlineSlope * reward_slopeLinearTrendline);

  return fullReward;
}

int fullSnapTest(const std::string& filename)
{
  // Create our factory:
  auto factory = JBrain::JBrainFactory::getInstance();
  factory->initialize(filename);

  // Create the gym runner:
  std::cout << "Creating GymSageRunner..." << std::endl;
  Experiment::GymSageRunner* sageRunner = Experiment::GymSageRunner::getInstance();
  // sageRunner->initialize(false, 4, 1, "CartPole-v1", "ppo_sage"); // useArgMax MUST be false:
  sageRunner->initialize(false, 8, 1, "LunarLander-v3", "ppo_sage_lunarLanderDiscrete");

  // Gather the experiment variables:
  YAML::Node expNode = factory->getExperimentConfig();
  std::string mainDataDirectory = expNode["MainDataDirectory"].as<std::string>();
  std::string expName = expNode["ExperimentName"].as<std::string>();
  unsigned int trainTrials = expNode["TrainTrials"].as<unsigned int>();
  unsigned int testTrials = expNode["TestTrials"].as<unsigned int>();
  unsigned int maxEpochs = expNode["MaximumEpochs"].as<unsigned int>();
  unsigned int startPopSize = expNode["StartingPopulationSize"].as<unsigned int>();
  bool useHDCMode = factory->getHDCModeSet();
  double maximumReward = expNode["MaximumReward"].as<double>();
  double popInfusionReward = expNode["PopulationInfusionReward"].as<double>();
  unsigned int popInfusionSize = expNode["PopulationInfusionSize"].as<unsigned int>();
  double reward_testNoOutput = expNode["Reward_TestNoOutput"].as<double>();
  double reward_testPercentGoodOutput = expNode["Reward_TestPercentGoodOutput"].as<double>();
  double reward_testPercentBadOutput = expNode["Reward_TestPercentBadOutput"].as<double>();
  double reward_testMaxScorePercent = expNode["Reward_TestMaxScorePercent"].as<double>();
  double reward_testBestMinusMin = expNode["Reward_TestBestMinusMin"].as<double>();
  double reward_testProcChurn = expNode["Reward_TestProcessingNeuronChurn"].as<double>();
  double reward_testInputChurn = expNode["Reward_TestInputNeuronChurn"].as<double>();
  double reward_fullTrendlineSlope = expNode["Reward_TrendlineScoreSlope"].as<double>();

  JBrain::JBrain_Snap* parent = nullptr;

  // Create our data directory. First: Extremely convoluted way to get current time:
  std::string dataDir;
  auto now = std::chrono::system_clock::now();
  time_t t = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(__unix__)
  localtime_r(&t, &tm);
#else
  localtime_s(&tm, &t);
#endif

  // Get the date-time as a string for a unique experiment directory name:
  char buffer[80];
  strftime(buffer, sizeof(buffer), "_%Y-%m-%d_%H-%M-%S", &tm);

  // Combine them into a single directory name:
  dataDir = mainDataDirectory + expName + std::string(buffer) + "\\";
  std::filesystem::create_directories(dataDir);

  // Copy our configuration file into the experiment directory:
  std::filesystem::copy_file(filename, dataDir + "experimentConfigFile.yaml");

  // Create the full experiment CSV:
  std::string csvFName = dataDir + "experimentSummary.csv";
  std::ofstream csvOut(csvFName.c_str(), std::ios_base::app);

  // Write the header and close. We re-open to append to prevent loss of data if the experiment crashes:
  csvOut << "epoch,bestScore,bestBrain,bestParent,trainingSeconds" << std::endl;
  csvOut.close();

  std::vector<JBrain::JBrain_Snap*> population;

  std::cout << "Getting " << startPopSize << " random sample brains..." << std::endl;
  for (unsigned int i = 0; i < startPopSize; ++i)
    population.push_back(factory->getRandomSnapBrain());

  std::vector<std::vector<double> >* m_obsData = nullptr;
  std::vector<std::vector<double> >* m_actData = nullptr;
  bool useCSVTrainingData = factory->getCSVTrainingDataProvided();
  if (useCSVTrainingData)
  {
    m_obsData = factory->getCSVObservations();
    m_actData = factory->getCSVActions();
  }
  
  // Test the full population:
  std::vector<double> allRewards;
  double maxReward;
  unsigned int maxIndex;
  for (unsigned int epoch = 0; epoch < maxEpochs; ++epoch)
  {
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Beginning training on epoch " << epoch + 1 << " / " << maxEpochs << std::endl;    
    allRewards.clear();

    for (auto brain : population)    
    {     
      // HDC Mode needs a different testing function:
      if (useHDCMode)
      {
        if (useCSVTrainingData)
          allRewards.push_back(testSnapJBrain_HDCCSV(brain, sageRunner, dataDir, testTrials, m_obsData, m_actData));
        else
          allRewards.push_back(testSnapJBrain_HDC(brain, sageRunner, dataDir, trainTrials, testTrials));
      }
      else
      {
        allRewards.push_back(testSnapJBrain(brain, sageRunner, dataDir, trainTrials, testTrials,
          reward_testNoOutput, reward_testPercentGoodOutput, reward_testPercentBadOutput,
          reward_testMaxScorePercent, reward_testBestMinusMin, reward_testProcChurn,
          reward_testInputChurn, reward_fullTrendlineSlope));
      }
      // std::cout << ".";
    }

    maxReward = -100.0;
    maxIndex = -1;
    for (unsigned int idx = 0; idx < allRewards.size(); ++idx)
    {
      if (allRewards[idx] > maxReward)
      {
        maxReward = allRewards[idx];
        maxIndex = idx;
      }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    // Gather the best:
    parent = population[maxIndex];
    std::cout << "\nMax reward: " << maxReward << " from brain " << parent->getName()
      << ". Full epoch took " << seconds << " seconds." << std::endl;

    // Write out to the csv file:
    csvOut.open(csvFName.c_str(), std::ios_base::app);
    csvOut << epoch << ","
      << maxReward << ","
      << parent->getName() << ","
      << parent->getParentName() << ","
      << seconds << std::endl;
    csvOut.close();

    // Delete all but the best:
    population.erase(population.begin() + maxIndex);
    for (auto brain : population)
      delete brain;

    // Create the new generation:
    population.clear();
    population = factory->getFullMutatedPopulation(parent);

    // Check if we're failing hard:
    if (maxReward < popInfusionReward)
    {
      std::cout << "\nPoor overall performance. " << maxReward << " < " << popInfusionReward
        << ". Adding " << popInfusionSize << " random individuals to the population." << std::endl;

      for (unsigned int i = 0; i < popInfusionSize; ++i)
        population.push_back(factory->getRandomSnapBrain());
    }

    // End the experiment early due to success:
    if (maxReward > maximumReward)
    {
      std::cout << "Reward greater than " << maximumReward << " achieved. Ending experiment." << std::endl;
      epoch = maxEpochs + 5;
    }
  }

  // Clean up:
  for (auto brain : population)
    delete brain;

  return 0;
}

int testTuples()
{
  typedef std::tuple<std::string, std::vector<std::string> > namePathTuple;

  std::vector<namePathTuple> testParams;

  testParams.push_back(namePathTuple("Hello", { "World.", "And stuff" }));
  testParams.push_back(namePathTuple("Goodbye", { "Friggen", "Whatever" }));

  std::cout << "testParams: " << std::endl;

  for (auto npt : testParams)
  {
    std::cout << std::get<0>(npt) << ": ";
    for (auto pathPart : std::get<1>(npt))
      std::cout << pathPart << " ";
    std::cout << std::endl;
  }

  return 0;
}

int main(int argc, char** argv)
{
  // return testObservationProcessor();
  // return testGymSageRunner_cartPole();
  // return testGymSageRunner_lunarLander();
  std::string configFileName = "snapBrainConfig.yaml";

  if (argc == 2)
    configFileName = std::string(argv[1]);
  // return testTuples();

  std::cout << "Beginning experiment described in file: " << configFileName << std::endl;
  return fullSnapTest(configFileName);
}
#endif
