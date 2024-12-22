#include "pch.h"
#include <cctype>
#include <string>
#include "Enums.h"

namespace CGP
{
	std::string DynamicProbabilityToString(DYNAMIC_PROBABILITY x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case DYNAMIC_PROBABILITY::ADD:
			retVal = "ADD";
			break;

		case DYNAMIC_PROBABILITY::SOLO:
			retVal = "SOLO";
			break;

		case DYNAMIC_PROBABILITY::MULTIPLY:
			retVal = "MULTIPLY";
			break;

		case DYNAMIC_PROBABILITY::UNUSED:
			retVal = "UNUSED";
			break;
		}

		return retVal;
	}

	DYNAMIC_PROBABILITY StringToDynamicProbability(std::string in)
	{
		DYNAMIC_PROBABILITY retVal = DYNAMIC_PROBABILITY::UNDEFINED;

		std::transform(in.begin(), in.end(), in.begin(), toupper);

		if (in == "ADD")
			retVal = DYNAMIC_PROBABILITY::ADD;
		else if (in == "SOLO")
			retVal = DYNAMIC_PROBABILITY::SOLO;
		else if (in == "MULTIPLY")
			retVal = DYNAMIC_PROBABILITY::MULTIPLY;
		else if (in == "UNUSED")
			retVal = DYNAMIC_PROBABILITY::UNUSED;

		return retVal;
	}

	std::string InputPreprocessingToString(INPUT_PREPROCESSING x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case INPUT_PREPROCESSING::NO_CHANGE:
			retVal = "NO_CHANGE";
			break;

		case INPUT_PREPROCESSING::BUCKETS:
			retVal = "BUCKETS";
			break;

		case INPUT_PREPROCESSING::NEGATIVE_VALUE_ADD:
			retVal = "NEGATIVE_VALUE_ADD";
			break;
		}
		return retVal;
	}

	INPUT_PREPROCESSING StringToInputPreprocessing(std::string in)
	{
		INPUT_PREPROCESSING retVal = INPUT_PREPROCESSING::UNDEFINED;

		std::transform(in.begin(), in.end(), in.begin(), toupper);

		if (in == "NO_CHANGE")
			retVal = INPUT_PREPROCESSING::NO_CHANGE;
		else if (in == "BUCKETS")
			retVal = INPUT_PREPROCESSING::BUCKETS;
		else if (in == "NEGATIVE_VALUE_ADD")
			retVal = INPUT_PREPROCESSING::NEGATIVE_VALUE_ADD;

		return retVal;
	}

	std::string ActivationFunctionToString(JNEURON_ACTIVATION_FUNCTION x)
	{
		std::string retVal = "NONE";
		switch (x)
		{
		case JNEURON_ACTIVATION_FUNCTION::TANH:
			retVal = "TANH";
			break;

		case JNEURON_ACTIVATION_FUNCTION::SIGMOID:
			retVal = "SIGMOID";
		}

		return retVal;
	}

	JNEURON_ACTIVATION_FUNCTION StringToActivationFunction(std::string in)
	{
		JNEURON_ACTIVATION_FUNCTION retVal = JNEURON_ACTIVATION_FUNCTION::NONE;

		std::transform(in.begin(), in.end(), in.begin(), toupper);

		if (in == "TANH")
			retVal = JNEURON_ACTIVATION_FUNCTION::TANH;
		else if (in == "SIGMOID")
			retVal = JNEURON_ACTIVATION_FUNCTION::SIGMOID;

		return retVal;
	}


	std::string UpdateEventToString(UPDATE_EVENT x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case UPDATE_EVENT::SCENARIO_RUN:
			retVal = "SCENARIO_RUN";
			break;

		case UPDATE_EVENT::BRAIN_OUTPUT:
			retVal = "BRAIN_OUTPUT";
			break;

		case UPDATE_EVENT::TIME_STEP:
			retVal = "TIME_STEP";
			break;
		}

		return retVal;
	}

	UPDATE_EVENT StringToUpdateEvent(std::string in)
	{
		UPDATE_EVENT retVal = UPDATE_EVENT::UNDEFINED;
				
		std::transform(in.begin(), in.end(), in.begin(), toupper);

		if (in == "SCENARIO_RUN")
			retVal = UPDATE_EVENT::SCENARIO_RUN;
		else if (in == "BRAIN_OUTPUT")
			retVal = UPDATE_EVENT::BRAIN_OUTPUT;
		else if (in == "TIME_STEP")
			retVal = UPDATE_EVENT::TIME_STEP;

		return retVal;
	}

	std::string JneuronSnapTypeToString(JNEURON_SNAP_TYPE x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case JNEURON_SNAP_TYPE::INPUT:
			retVal = "INPUT";
			break;
		case JNEURON_SNAP_TYPE::PROCESSING:
			retVal = "PROCESSING";
			break;
		case JNEURON_SNAP_TYPE::OUTPUT:
			retVal = "OUTPUT";
			break;
		}

		return retVal;
	}

	JNEURON_SNAP_TYPE StringToJneuronSnapType(std::string in)
	{
		JNEURON_SNAP_TYPE retVal = JNEURON_SNAP_TYPE::UNDEFINED;

		std::transform(in.begin(), in.end(), in.begin(), toupper);

		if (in == "INPUT")
			retVal = JNEURON_SNAP_TYPE::INPUT;
		else if (in == "PROCESSING")
			retVal = JNEURON_SNAP_TYPE::PROCESSING;
		else if (in == "OUTPUT")
			retVal = JNEURON_SNAP_TYPE::OUTPUT;

		return retVal;
	}

	std::string CGPInputToString(CGP_INPUT x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case CGP_INPUT::SAGE_MATCH_PERCENT:
			retVal = "SAGE_MATCH_PERCENT";
			break;
		case CGP_INPUT::CURRENT_WEIGHT:
			retVal = "CURRENT_WEIGHT";
			break;
		case CGP_INPUT::STRONGEST_INPUT_XYZ:
			retVal = "STRONGEST_INPUT_XYZ";
			break;
		case CGP_INPUT::STRONGEST_INPUT_DISTANCE:
			retVal = "STRONGEST_INPUT_DISTANCE";
			break;
		case CGP_INPUT::STRONGEST_INPUT_IS_OBSERVATION_AXON:
			retVal = "STRONGEST_INPUT_IS_OBSERVATION_AXON";
			break;
		case CGP_INPUT::STRONGEST_INPUT_VALUE:
			retVal = "STRONGEST_INPUT_VALUE";
			break;
		case CGP_INPUT::NEAREST_AXON_XYZ:
			retVal = "NEAREST_AXON_XYZ";
			break;
		case CGP_INPUT::NEAREST_AXON_DISTANCE:
			retVal = "NEAREST_AXON_DISTANCE";
			break;
		case CGP_INPUT::NEAREST_AXON_IS_OBSERVATION_AXON:
			retVal = "NEAREST_AXON_IS_OBSERVATION_AXON";
			break;
		case CGP_INPUT::NEAREST_AXON_IS_PART_OF_SAME_NEURON:
			retVal = "NEAREST_AXON_IS_PART_OF_SAME_NEURON";
			break;
		case CGP_INPUT::NEAREST_DENDRITE_XYZ:
			retVal = "NEAREST_DENDRITE_XYZ";
			break;
		case CGP_INPUT::NEAREST_DENDRITE_IS_PART_OF_SAME_NEURON:
			retVal = "NEAREST_DENDRITE_IS_PART_OF_SAME_NEURON";
			break;
		case CGP_INPUT::DENDRITE_TYPE:
			retVal = "DENDRITE_TYPE";
			break;
		case CGP_INPUT::INPUT_MAGNITUDE:
			retVal = "INPUT_MAGNITUDE";
			break;
		case CGP_INPUT::CURRENT_LENGTH:
			retVal = "CURRENT_LENGTH";
			break;
		case CGP_INPUT::NEURON_AGE:
			retVal = "NEURON_AGE";
			break;
		case CGP_INPUT::NEURON_HEALTH:
			retVal = "NEURON_HEALTH";
			break;
		case CGP_INPUT::PERCENTAGE_FIRE:
			retVal = "PERCENTAGE_FIRE";
			break;
		case CGP_INPUT::PERCENTAGE_BRAIN_FIRE:
			retVal = "PERCENTAGE_BRAIN_FIRE";
			break;
		case CGP_INPUT::EXPECTED_OUTPUT_DIFF:
			retVal = "EXPECTED_OUTPUT_DIFF";
			break;
		}

		return retVal;
	}

	CGP_INPUT StringToCGPInput(std::string in)
	{
		CGP_INPUT retVal = CGP_INPUT::UNDEFINED;
		std::transform(in.begin(), in.end(), in.begin(), toupper);

		if (in == "SAGE_MATCH_PERCENT")
			retVal = CGP_INPUT::SAGE_MATCH_PERCENT;
		else if (in == "CURRENT_WEIGHT")
			retVal = CGP_INPUT::CURRENT_WEIGHT;
		else if (in == "STRONGEST_INPUT_XYZ")
			retVal = CGP_INPUT::STRONGEST_INPUT_XYZ;
		else if (in == "STRONGEST_INPUT_DISTANCE")
			retVal = CGP_INPUT::STRONGEST_INPUT_DISTANCE;
		else if (in == "STRONGEST_INPUT_IS_OBSERVATION_AXON")
			retVal = CGP_INPUT::STRONGEST_INPUT_IS_OBSERVATION_AXON;
		else if (in == "STRONGEST_INPUT_VALUE")
			retVal = CGP_INPUT::STRONGEST_INPUT_VALUE;
		else if (in == "NEAREST_AXON_XYZ")
			retVal = CGP_INPUT::NEAREST_AXON_XYZ;
		else if (in == "NEAREST_AXON_DISTANCE")
			retVal = CGP_INPUT::NEAREST_AXON_DISTANCE;
		else if (in == "NEAREST_AXON_IS_OBSERVATION_AXON")
			retVal = CGP_INPUT::NEAREST_AXON_IS_OBSERVATION_AXON;
		else if (in == "NEAREST_AXON_IS_PART_OF_SAME_NEURON")
			retVal = CGP_INPUT::NEAREST_AXON_IS_PART_OF_SAME_NEURON;
		else if (in == "NEAREST_DENDRITE_XYZ")
			retVal = CGP_INPUT::NEAREST_DENDRITE_XYZ;
		else if (in == "NEAREST_DENDRITE_IS_PART_OF_SAME_NEURON")
			retVal = CGP_INPUT::NEAREST_DENDRITE_IS_PART_OF_SAME_NEURON;
		else if (in == "DENDRITE_TYPE")
			retVal = CGP_INPUT::DENDRITE_TYPE;
		else if (in == "INPUT_MAGNITUDE")
			retVal = CGP_INPUT::INPUT_MAGNITUDE;
		else if (in == "CURRENT_LENGTH")
			retVal = CGP_INPUT::CURRENT_LENGTH;
		else if (in == "NEURON_AGE")
			retVal = CGP_INPUT::NEURON_AGE;
		else if (in == "NEURON_HEALTH")
			retVal = CGP_INPUT::NEURON_HEALTH;
		else if (in == "PERCENTAGE_FIRE")
			retVal = CGP_INPUT::PERCENTAGE_FIRE;
		else if (in == "PERCENTAGE_BRAIN_FIRE")
			retVal = CGP_INPUT::PERCENTAGE_BRAIN_FIRE;
		else if (in == "EXPECTED_OUTPUT")
			retVal = CGP_INPUT::EXPECTED_OUTPUT_DIFF;

		return retVal;
	}

	std::string CGPOutputToString(CGP_OUTPUT x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case CGP_OUTPUT::LOCATION:
			retVal = "LOCATION";
			break;
		case CGP_OUTPUT::STRONGEST_INPUT_CLOSER_FURTHER:
			retVal = "STRONGEST_INPUT_CLOSER_FURTHER";
			break;
		case CGP_OUTPUT::NEAREST_AXON_CLOSER_FURTHER:
			retVal = "NEAREST_AXON_CLOSER_FURTHER";
			break;
		case CGP_OUTPUT::RANDOM_MOVEMENT_THRESHOLD:
			retVal = "RANDOM_MOVEMENT_THRESHOLD";
		case CGP_OUTPUT::NEAREST_DENDRITE_CLOSER_FURTHER:
			retVal = "NEAREST_DENDRITE_CLOSER_FURTHER";
			break;
		case CGP_OUTPUT::CLOSER_TO_STRONGEST_INPUT:
			retVal = "CLOSER_TO_STRONGEST_INPUT";
			break;
		case CGP_OUTPUT::CLOSER_TO_NEAREST_AXON:
			retVal = "CLOSER_TO_NEAREST_AXON";
			break;
		case CGP_OUTPUT::HEALTH:
			retVal = "HEALTH";
			break;
		case CGP_OUTPUT::HEALTH_INCREASE_DECREASE:
			retVal = "HEALTH_INCREASE_DECREASE";
			break;
		case CGP_OUTPUT::WEIGHT:
			retVal = "WEIGHT";
			break;
		case CGP_OUTPUT::WEIGHT_HIGHER_LOWER:
			retVal = "WEIGHT_HIGHER_LOWER";
			break;
		}

		return retVal;
	}
	
	CGP_OUTPUT StringToCGPOutput(std::string in)
	{
		CGP_OUTPUT retVal = CGP_OUTPUT::UNDEFINED;
		std::transform(in.begin(), in.end(), in.begin(), toupper);

		if (in == "LOCATION")
			retVal = CGP_OUTPUT::LOCATION;
		else if (in == "STRONGEST_INPUT_CLOSER_FURTHER")
			retVal = CGP_OUTPUT::STRONGEST_INPUT_CLOSER_FURTHER;
		else if (in == "NEAREST_AXON_CLOSER_FURTHER")
			retVal = CGP_OUTPUT::NEAREST_AXON_CLOSER_FURTHER;
		else if (in == "RANDOM_MOVEMENT_THRESHOLD")
			retVal = CGP_OUTPUT::RANDOM_MOVEMENT_THRESHOLD;
		else if (in == "NEAREST_DENDRITE_CLOSER_FURTHER")
			retVal = CGP_OUTPUT::NEAREST_DENDRITE_CLOSER_FURTHER;
		else if (in == "CLOSER_TO_STRONGEST_INPUT")
			retVal = CGP_OUTPUT::CLOSER_TO_STRONGEST_INPUT;
		else if (in == "CLOSER_TO_NEAREST_AXON")
			retVal = CGP_OUTPUT::CLOSER_TO_NEAREST_AXON;
		else if (in == "HEALTH")
			retVal = CGP_OUTPUT::HEALTH;
		else if (in == "HEALTH_INCREASE_DECREASE")
			retVal = CGP_OUTPUT::HEALTH_INCREASE_DECREASE;
		else if (in == "WEIGHT")
			retVal = CGP_OUTPUT::WEIGHT;
		else if (in == "WEIGHT_HIGHER_LOWER")
			retVal = CGP_OUTPUT::WEIGHT_HIGHER_LOWER;

		return retVal;
	}

	std::string FitnessCollapseFunctionToString(FITNESS_COLLAPSE_FUNCTION x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case FITNESS_COLLAPSE_FUNCTION::MEAN:
			retVal = "MEAN";
			break;
		case FITNESS_COLLAPSE_FUNCTION::MEDIAN:
			retVal = "MEDIAN";
			break;
		case FITNESS_COLLAPSE_FUNCTION::AVG_OF_MEAN_AND_MEDIAN:
			retVal = "AVG_OF_MEAN_AND_MEDIAN";
			break;
		case FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN:
			retVal = "MIN_OF_MEAN_AND_MEDIAN";
			break;
		case FITNESS_COLLAPSE_FUNCTION::MINIMUM:
			retVal = "MINIMUM";
			break;
		case FITNESS_COLLAPSE_FUNCTION::MAXIMUM:
			retVal = "MAXIMUM";
			break;
		case FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN_PLUS_MIN:
			retVal = "MIN_OF_MEAN_AND_MEDIAN_PLUS_MIN";
			break;
		case FITNESS_COLLAPSE_FUNCTION::UNDEFINED:
			retVal = "UNDEFINED";
			break;
		}

		return retVal;
	}

	FITNESS_COLLAPSE_FUNCTION StringToFitnessCollapseFunction(std::string in)
	{
		std::transform(in.begin(), in.end(), in.begin(), toupper);
		FITNESS_COLLAPSE_FUNCTION retVal = FITNESS_COLLAPSE_FUNCTION::UNDEFINED;
		
		if (in == "MEAN")
			retVal = FITNESS_COLLAPSE_FUNCTION::MEAN;
		else if (in == "MEDIAN")
			retVal = FITNESS_COLLAPSE_FUNCTION::MEDIAN;
		else if (in == "AVG_OF_MEAN_AND_MEDIAN")
			retVal = FITNESS_COLLAPSE_FUNCTION::AVG_OF_MEAN_AND_MEDIAN;
		else if (in == "MIN_OF_MEAN_AND_MEDIAN")
			retVal = FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN;
		else if (in == "MINIMUM")
			retVal = FITNESS_COLLAPSE_FUNCTION::MINIMUM;
		else if (in == "MAXIMUM")
			retVal = FITNESS_COLLAPSE_FUNCTION::MAXIMUM;
		else if (in == "MIN_OF_MEAN_AND_MEDIAN_PLUS_MIN")
			retVal = FITNESS_COLLAPSE_FUNCTION::MIN_OF_MEAN_AND_MEDIAN_PLUS_MIN;

		return retVal;
	}

	std::string FunctionTypeToString(FUNCTION_TYPE x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case FUNCTION_TYPE::BOOL3_BOOL:
			retVal = "BOOL3_BOOL";
			break;
		case FUNCTION_TYPE::DOUBLE3_DOUBLE:
			retVal = "DOUBLE3_DOUBLE";
			break;
		case FUNCTION_TYPE::DOUBLE_DOUBLE:
			retVal = "DOUBLE_DOUBLE";
			break;
		case FUNCTION_TYPE::UNDEFINED:
			retVal = "UNDEFINED";
			break;
		}

		return retVal;
	}

	FUNCTION_TYPE StringToFunctionType(std::string in)
	{
		std::transform(in.begin(), in.end(), in.begin(), toupper);
		FUNCTION_TYPE retVal = FUNCTION_TYPE::UNDEFINED;

		if (in == "BOOL3_BOOL")
			retVal = FUNCTION_TYPE::BOOL3_BOOL;
		else if (in == "DOUBLE3_DOUBLE")
			retVal = FUNCTION_TYPE::DOUBLE3_DOUBLE;
		else if (in == "DOUBLE_DOUBLE")
			retVal = FUNCTION_TYPE::DOUBLE_DOUBLE;

		return retVal;
	}

	std::string GeneTypeToString(GENE_TYPE x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case GENE_TYPE::INPUT:
			retVal = "INPUT";
			break;
		case GENE_TYPE::PROCESSING:
			retVal = "PROCESSING";
			break;
		case GENE_TYPE::OUTPUT:
			retVal = "OUTPUT";
			break;
		case GENE_TYPE::UNDEFINED:
			retVal = "UNDEFINED";
			break;
		}

		return retVal;
	}

	GENE_TYPE StringToGeneType(std::string in)
	{
		std::transform(in.begin(), in.end(), in.begin(), toupper);
		GENE_TYPE retVal = GENE_TYPE::UNDEFINED;

		if (in == "INPUT")
			retVal = GENE_TYPE::INPUT;
		else if (in == "PROCESSING")
			retVal = GENE_TYPE::PROCESSING;
		else if (in == "OUTPUT")
			retVal = GENE_TYPE::OUTPUT;

		return retVal;
	}

	std::string CGPTypeToString(CGP_TYPE x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case CGP_TYPE::GENERAL:
			retVal = "GENERAL";
			break;
		case CGP_TYPE::CGPFFANN:
			retVal = "CGPFFANN";
			break;
		case CGP_TYPE::UNDEFINED:
			retVal = "UNDEFINED";
			break;
		}

		return retVal;
	}
	CGP_TYPE StringToCGPType(std::string in)
	{
		std::transform(in.begin(), in.end(), in.begin(), toupper);
		CGP_TYPE retVal = CGP_TYPE::UNDEFINED;

		if (in == "GENERAL")
			retVal = CGP_TYPE::GENERAL;
		else if (in == "CGPFFANN")
			retVal = CGP_TYPE::CGPFFANN;

		return retVal;
	}

	std::string MutationStrategyToString(MUTATION_STRATEGY x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case MUTATION_STRATEGY::PROBABILISTIC_PERGENE:
			retVal = "PROBABILISTIC_PERGENE";
			break;
		case MUTATION_STRATEGY::PROBABILISTIC_PERVALUE:
			retVal = "PROBABILISTIC_PERVALUE";
			break;
		case MUTATION_STRATEGY::ACTIVE_GENE:
			retVal = "ACTIVE_GENE";
			break;
		case MUTATION_STRATEGY::UNDEFINED:
			retVal = "UNDEFINED";
			break;
		}

		return retVal;
	}

	MUTATION_STRATEGY StringToMutationStrategy(std::string in)
	{
		std::transform(in.begin(), in.end(), in.begin(), toupper);
		MUTATION_STRATEGY retVal = MUTATION_STRATEGY::UNDEFINED;

		if (in == "PROBABILISTIC_PERGENE")
			retVal = MUTATION_STRATEGY::PROBABILISTIC_PERGENE;
		else if (in == "PROBABILISTIC_PERVALUE")
			retVal = MUTATION_STRATEGY::PROBABILISTIC_PERVALUE;
		else if (in == "ACTIVE_GENE")
			retVal = MUTATION_STRATEGY::ACTIVE_GENE;

		return retVal;
	}

	std::string GymTypeToString(GYM_TYPE x)
	{
		std::string retVal = "UNDEFINED";
		switch (x)
		{
		case GYM_TYPE::PyGYM:
			retVal = "PYGYM";
			break;
		case GYM_TYPE::CGYM_1D:
			retVal = "CGYM_1D";
			break;
		case GYM_TYPE::CGYM_2D:
			retVal = "CGYM_2D";
			break;
		case GYM_TYPE::UNDEFINED:
			retVal = "UNDEFINED";
			break;
		}

		return retVal;
	}
	GYM_TYPE StringToGymType(std::string in)
	{
		std::transform(in.begin(), in.end(), in.begin(), toupper);
		GYM_TYPE retVal = GYM_TYPE::UNDEFINED;

		if (in == "PYGYM")
			retVal = GYM_TYPE::PyGYM;
		else if (in == "CGYM_1D")
			retVal = GYM_TYPE::CGYM_1D;
		else if (in == "CGYM_2D")
			retVal = GYM_TYPE::CGYM_2D;

		return retVal;
	}

}