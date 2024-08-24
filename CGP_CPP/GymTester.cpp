#include "pch.h"
#include "GymTester.h"
#include <algorithm>
#include <numeric>

// Debug:
#include <iostream>

extern unsigned int DEBUG_LEVEL;

namespace CGP
{
    GymTester* GymTester::getInstance()
    {
        // Not built to deal with the potential of multiple threads:
        static GymTester* instance = new GymTester();

        return instance;
    }

    GymTester::GymTester()
        : AbstractTesterClass(),
        m_gymModule(nullptr),
        m_env(nullptr),
        m_envReset(nullptr),
        m_envStep(nullptr),
        m_envClose(nullptr),
        m_envRender(nullptr),
        m_gcModule(nullptr),
        m_gcCollectFunction(nullptr),
        m_observationSize(-1),
        m_actionSize(-1),
        m_collapseFunction(FITNESS_COLLAPSE_FUNCTION::UNDEFINED),
        m_useArgMax(false),
        m_initialized(false)         
    {          
    }
    
    void GymTester::initializePython()
	{
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		Py_Initialize();

        // 2 lines to load our custom environment:
        PyRun_SimpleString("from gym.envs.registration import register");
        PyRun_SimpleString("register(id='DungeonRoomEnv', entry_point='gym.envs.customEnvs.DungeonRoomEnv:DungeonRoomEnv', max_episode_steps=300, reward_threshold=105.0)");

		// Get the gym module:        
        PyObject* pGymName = PyUnicode_FromString("gym");
        m_gymModule = PyImport_Import(pGymName);

        if (m_gymModule == nullptr)
        {
            PyErr_Print();
            exit(-1);
        }

        // Get the make function:
        PyObject* pMakeFunc = PyObject_GetAttrString(m_gymModule, "make");

        if (pMakeFunc == nullptr)
        {
            PyErr_Print();
            exit(-1);
        }

        PyObject* pMakeArgs = PyTuple_New(1);
        PyTuple_SetItem(pMakeArgs, 0, PyUnicode_FromString(m_envName.c_str()));

        m_env = PyObject_CallObject(pMakeFunc, pMakeArgs);

        // Get the functions we'll need later:
        m_envReset = PyObject_GetAttrString(m_env, "reset");
        
        if (m_envReset == nullptr)
        {
            std::cout << "reset error: ";
            PyErr_Print();
            exit(-1);
        }
        
        m_envStep = PyObject_GetAttrString(m_env, "step");
        
        if (m_envStep == nullptr)
        {
            std::cout << "step error: ";
            PyErr_Print();
            exit(-1);
        }
        
        m_envClose = PyObject_GetAttrString(m_env, "close");

        if (m_envClose == nullptr)
        {
            std::cout << "close error: ";
            PyErr_Print();
            exit(-1);
        }
        m_envRender = PyObject_GetAttrString(m_env, "render");

        if (m_envRender == nullptr)
        {
            std::cout << "render error: ";
            PyErr_Print();
            exit(-1);
        }

        Py_CLEAR(pGymName);
        Py_CLEAR(pMakeFunc);
        Py_CLEAR(pMakeArgs);

        // Get the garbage collection module and function:
        PyObject* pGCName = PyUnicode_FromString("gc");
        m_gcModule = PyImport_Import(pGCName);
        m_gcCollectFunction = PyObject_GetAttrString(m_gcModule, "collect");

        Py_CLEAR(pGCName);
    }

	void GymTester::finalizePython()
	{
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Clear all Gym pointers:
        Py_CLEAR(m_envReset);
        Py_CLEAR(m_envStep);
        Py_CLEAR(m_envClose);
        Py_CLEAR(m_envRender);

        Py_CLEAR(m_env);

        // Clear garbage collection pointers:
        Py_CLEAR(m_gcCollectFunction);
        Py_CLEAR(m_gcModule);
        
        Py_Finalize();
	}

	void GymTester::initialize(
        std::string envName, FITNESS_COLLAPSE_FUNCTION fitFunc,
		bool useArgMax, int observationSize, int actionSize)		
	{
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Only go through initialization if it is the first time:
        if (!m_initialized)
        {
            m_envName = envName;
            m_collapseFunction = fitFunc;
            m_useArgMax = useArgMax;
            m_observationSize = observationSize;
            m_actionSize = actionSize;
            m_initialized = true;

            initializePython();
        }
	}

	GymTester::~GymTester()
	{   
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

		finalizePython();
	}

    double GymTester::renderIndividual(CGP::AbstractCGPIndividual* ind, int msRenderDelay,
                                       unsigned int maxSteps)
    {        
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Get a score while requesting a render:
        double retScore = getIndScore(ind, msRenderDelay, 1, maxSteps);
             
        return retScore;
    }

	double GymTester::getIndScore(AbstractCGPIndividual* ind, int msRenderDelay,
                                  int runsPerFitness, unsigned int maxStepsPerRun)
	{
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        bool done;
        double fullReward;
        std::vector<double> allRewards;
        std::vector<double> observation;
        std::vector<double> action;
        int maxElemIdx = -1;
        unsigned int stepNumber = 0;
        int actionSize = m_actionSize;
        if (m_useArgMax)
            actionSize = 1;
        
        for (int i = 0; i < runsPerFitness; ++i)
        {
            done = false;
            fullReward = 0.0;
            stepNumber = 0;
            
            // Reset the environment and get the initial observation:
            PyObject* pInitialObsArray = PyObject_CallNoArgs(m_envReset);            
            
            // The Gym library seems to be changing/updating and starting to return
            // a tuple from reset. We need the array at the front:
            PyObject* pInitialObsListFunc = nullptr;
            PyObject* pActualObsArray = nullptr;

            if (PyTuple_Check(pInitialObsArray))
            {
                pActualObsArray = PyTuple_GetItem(pInitialObsArray, 0);
                pInitialObsListFunc = PyObject_GetAttrString(pActualObsArray, "tolist");
            }
            else
            {
                pInitialObsListFunc = PyObject_GetAttrString(pInitialObsArray, "tolist");
            }

            if (pInitialObsListFunc == nullptr)
            {
                PyErr_Print();
                exit(-1);
            }
            
            PyObject* pInitialObsList = PyObject_CallNoArgs(pInitialObsListFunc);

            observation.clear();
            for (int j = 0; j < m_observationSize; ++j)
            {
                // PyList_GetItem returns don't need to be decref'd; it returns a borrowed ref:
                observation.push_back(PyFloat_AsDouble(PyList_GetItem(pInitialObsList, j)));
            }

            while (!done)
            {
                if (msRenderDelay >= 0)
                {
                    PyObject_CallNoArgs(m_envRender);                    
                    std::this_thread::sleep_for(std::chrono::milliseconds(msRenderDelay));
                }

                // Have the individual calculate its outputs:
                action = ind->calculateOutputs(observation);                

                // If we're not using argMax, we need the action as a PyList:
                PyObject* pActionList = nullptr;
                if (!m_useArgMax)
                {
                    pActionList = PyList_New(actionSize);

                    for (int j = 0; j < actionSize; ++j)
                    {
                        PyList_SetItem(pActionList, j, PyFloat_FromDouble(action[j]));
                    }
                }

                // Put the action in a tuple for passing to the step function:
                PyObject* pActionTuple = PyTuple_New(1);
                
                if (m_useArgMax)
                {
                    maxElemIdx = static_cast<int>(std::max_element(action.begin(),
                        action.end()) - action.begin());
                    PyTuple_SetItem(pActionTuple, 0, PyLong_FromLong(maxElemIdx));
                }
                else
                {
                    PyTuple_SetItem(pActionTuple, 0, pActionList);
                }

                // Call the step function:
                PyObject* pStepReturn = PyObject_CallObject(m_envStep, pActionTuple);
                if (pStepReturn == nullptr)
                {
                    std::cout << "pStepReturn Null." << std::endl;
                    PyErr_Print();
                    assert(false);
                }

                // Get the 3 values we need from the step return:
                // 1. Reward:
                PyObject* pRewardReturn = PyTuple_GetItem(pStepReturn, 1);
                fullReward += PyFloat_AsDouble(pRewardReturn);

                // 2. Done:
                PyObject* pDoneReturn = PyTuple_GetItem(pStepReturn, 2);
                done = (Py_True == pDoneReturn) || (stepNumber > maxStepsPerRun);

                // 3. Next Observation:
                PyObject* pObsArray = PyTuple_GetItem(pStepReturn, 0);
                PyObject* pObsArrayListFunc = PyObject_GetAttrString(pObsArray, "tolist");
                PyObject* pObsList = PyObject_CallNoArgs(pObsArrayListFunc);

                observation.clear();
                for (int j = 0; j < m_observationSize; ++j)  // 8 is env obs-size
                {
                    // PyList_GetItem returns don't need to be decref'd; it is a borrowed ref:
                    observation.push_back(PyFloat_AsDouble(PyList_GetItem(pObsList, j)));
                }

                Py_CLEAR(pObsList);
                Py_CLEAR(pObsArrayListFunc);
                Py_CLEAR(pStepReturn);
                Py_CLEAR(pActionTuple);
                ++stepNumber;
            }
            
            allRewards.push_back(fullReward);
            
            Py_CLEAR(pInitialObsList);
            Py_CLEAR(pInitialObsListFunc);
            Py_CLEAR(pInitialObsArray);
        }

        return applyCollapseFunction(allRewards);
	}

	std::vector<double> GymTester::getPopulationScores(
		std::vector<CGP::AbstractCGPIndividual* > population,
        int runsPerFitness, unsigned int maxStepsPerRun)
	{
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        std::vector<double> retVal;        
        for (auto ind : population)
        {
            retVal.push_back(getIndScore(ind, -1, runsPerFitness, maxStepsPerRun));
        }
        
        // Manually invoke the garbage collection:
        PyObject* pGCReturn = PyObject_CallNoArgs(m_gcCollectFunction);
        auto objectsCollected = PyLong_AsLong(pGCReturn);
        if (objectsCollected > 0)
        {
            //std::cout << "***********\nObjects garbage collected: " << objectsCollected
            //    << "\n************" << std::endl;
        }

        Py_CLEAR(pGCReturn);
        
        return retVal;
    }

    double GymTester::applyCollapseFunction(std::vector<double> indRewards)
    {
        if (DEBUG_LEVEL > 4)
            std::cout << "At top of \"" << classname() << ":" << __func__ << "\"" << std::endl;

        // Only one value? Return it:
        if (indRewards.size() == 1)
        {
            return indRewards[0];
        }

        // For now, just go with min of median and mean:
        // Sort to get median:
        std::sort(indRewards.begin(), indRewards.end());
        
        double median = indRewards[indRewards.size() / 2];
        if (indRewards.size() % 2 == 0)
        {
            // If the length is even, need to include the value before:
            median = (median + indRewards[(indRewards.size() / 2) - 1]) / 2.0;
        }

        double mean = std::accumulate(indRewards.begin(), indRewards.end(), 0.0) / double(indRewards.size());

        // Return the minimum of the two:
        double retVal = median;
        if (mean < median)
            retVal = mean;

        return retVal;
    }
}; // End namespace CGP
