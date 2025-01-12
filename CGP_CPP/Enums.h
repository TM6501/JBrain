#pragma once
#include <string>
#include <algorithm>

// StringToXXX functions don't take const refs because then we just
// need to make a copy, anyway.

namespace CGP
{
	// If/How HDC learns:
	enum class HDC_LEARN_MODE : int {
		FULL, // Always update our understanding
		WRONG, // Update if we answered incorrectly
		NONE,  // No updates (probably testing)
		UNDEFINED
	};
	
	std::string HDCLearnModeToString(HDC_LEARN_MODE x);
	HDC_LEARN_MODE StringToHDCLearnMode(std::string in);

	// The dynamic overall probability type:
	enum class DYNAMIC_PROBABILITY : int {
		ADD,
		SOLO,
		MULTIPLY,
		UNUSED,
		UNDEFINED
	};

	std::string DynamicProbabilityToString(DYNAMIC_PROBABILITY x);
	DYNAMIC_PROBABILITY StringToDynamicProbability(std::string in);

	// The input preprocessing type:
	enum class INPUT_PREPROCESSING : int {
		UNDEFINED,
		NO_CHANGE,
		BUCKETS,
		NEGATIVE_VALUE_ADD
	};

	std::string InputPreprocessingToString(INPUT_PREPROCESSING x);
	INPUT_PREPROCESSING StringToInputPreprocessing(std::string in);

	// The activation functions that can be applied to neuron outputs:
	enum class JNEURON_ACTIVATION_FUNCTION :int {
		NONE,
		TANH,
		SIGMOID
	};

	std::string ActivationFunctionToString(JNEURON_ACTIVATION_FUNCTION x);
	JNEURON_ACTIVATION_FUNCTION StringToActivationFunction(std::string in);

	// Update frequency determines the triggering event that causes
	// the update programs to run and change the JBrain's components:
	enum class UPDATE_EVENT:int {
		SCENARIO_RUN,
		BRAIN_OUTPUT,
		TIME_STEP,
		UNDEFINED
	};

	std::string UpdateEventToString(UPDATE_EVENT x);
	UPDATE_EVENT StringToUpdateEvent(std::string in);	
	
	// The possible types of neurons used by the snap paradigm:
	enum class JNEURON_SNAP_TYPE : int {
		UNDEFINED,
		INPUT,
		PROCESSING,
		OUTPUT
	};

	std::string JneuronSnapTypeToString(JNEURON_SNAP_TYPE x);
	JNEURON_SNAP_TYPE StringToJneuronSnapType(std::string in);

	// The inputs that are used by the CGP programs that update
	// the brain's characteristics:
	enum class CGP_INPUT:int {
		SAGE_MATCH_PERCENT,
		CURRENT_WEIGHT,
		STRONGEST_INPUT_XYZ,
		STRONGEST_INPUT_DISTANCE,
		STRONGEST_INPUT_IS_OBSERVATION_AXON,
		STRONGEST_INPUT_VALUE,
		NEAREST_AXON_XYZ,
		NEAREST_AXON_DISTANCE,
		NEAREST_AXON_IS_OBSERVATION_AXON,
		NEAREST_AXON_IS_PART_OF_SAME_NEURON,
		NEAREST_DENDRITE_XYZ,
		NEAREST_DENDRITE_IS_PART_OF_SAME_NEURON,
		DENDRITE_TYPE,
		INPUT_MAGNITUDE,
		CURRENT_LENGTH,
		NEURON_AGE,
		NEURON_HEALTH,
		PERCENTAGE_FIRE,
		PERCENTAGE_BRAIN_FIRE,
		EXPECTED_OUTPUT_DIFF,
		UNDEFINED
	};

	std::string CGPInputToString(CGP_INPUT x);
	CGP_INPUT StringToCGPInput(std::string in);

	// The outputs that are used by the CGP programs that update the
	// brain's characteristics:
	enum class CGP_OUTPUT:int {
		LOCATION,
		STRONGEST_INPUT_CLOSER_FURTHER,
		NEAREST_AXON_CLOSER_FURTHER,
		RANDOM_MOVEMENT_THRESHOLD,
		NEAREST_DENDRITE_CLOSER_FURTHER,
		CLOSER_TO_STRONGEST_INPUT,
		CLOSER_TO_NEAREST_AXON,
		HEALTH,
		HEALTH_INCREASE_DECREASE,
		WEIGHT,
		WEIGHT_HIGHER_LOWER,
		UNDEFINED
	};

	std::string CGPOutputToString(CGP_OUTPUT x);
	CGP_OUTPUT StringToCGPOutput(std::string in);

	// OpenAI Gym environments include randomized starting positions.
	// Often times, we want to run a single individual through multiple
	// times to get a sense of its actual capabilities and minimize the
	// effect of luck. These functions define how we'll take a series of
	// individual fitness values and collapse them to a single fitness that
	// can be compared with other individuals:
	enum class FITNESS_COLLAPSE_FUNCTION:int {
		MEAN,
		MEDIAN,
		AVG_OF_MEAN_AND_MEDIAN,
		MIN_OF_MEAN_AND_MEDIAN,
		MINIMUM,
		MAXIMUM,
		// Need to be careful with this one. It will help push evolution in
		// the right direction, but it changes the maximum possible fitness.
		// This needs to be accounted for in the fitness values:
		MIN_OF_MEAN_AND_MEDIAN_PLUS_MIN,
		UNDEFINED
	};

	std::string FitnessCollapseFunctionToString(FITNESS_COLLAPSE_FUNCTION x);
	FITNESS_COLLAPSE_FUNCTION StringToFitnessCollapseFunction(std::string in);

	// These are the types of functions the CGP uses in its neurons.
	// Individuals will need to hold this information so that they can
	// read and write themselves to/from file:
	enum class FUNCTION_TYPE:int {
		BOOL3_BOOL,  // 3 Bool in, 1 bool out.
		DOUBLE3_DOUBLE,  // 3 double in, double out.
		DOUBLE_DOUBLE, // Neuron activation. Double in, double out.
		UNDEFINED
	};

	std::string FunctionTypeToString(FUNCTION_TYPE x);
	FUNCTION_TYPE StringToFunctionType(std::string in);

	enum class GENE_TYPE:int {
		INPUT,
		PROCESSING,
		OUTPUT,
		UNDEFINED
	};

	std::string GeneTypeToString(GENE_TYPE x);
	GENE_TYPE StringToGeneType(std::string in);

	enum class CGP_TYPE:int {
		GENERAL,
		CGPFFANN,
		UNDEFINED
	};
	
	std::string CGPTypeToString(CGP_TYPE x);
	CGP_TYPE StringToCGPType(std::string in);

	enum class MUTATION_STRATEGY:int {
		PROBABILISTIC_PERGENE,  // Mutate a certain percentage of genes.
		PROBABILISTIC_PERVALUE, // Mutate a certain percentage of values 
								// in the genetic code.
		ACTIVE_GENE,  // Mutate until a certain number of active genes are mutated.
		UNDEFINED
	};

	std::string MutationStrategyToString(MUTATION_STRATEGY x);
	MUTATION_STRATEGY StringToMutationStrategy(std::string in);	

	enum class GYM_TYPE :int {
		PyGYM,  // The traditional Gym worlds created in Python.
		// The RL worlds created in C++ which adhere to the Python Gym interface.
		// 1D / 2D refers to the shape of the observations that the environment provides.
		CGYM_1D,
		CGYM_2D,
		UNDEFINED
	};

	std::string GymTypeToString(GYM_TYPE x);
	GYM_TYPE StringToGymType(std::string in);
}