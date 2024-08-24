#pragma once
#include <vector>
#include <iostream>

// JBrain components will all be structures. They are meant to be used
// and manipulated freely by the JBrains that create them.
namespace JBrain
{
	struct JBrainComponent
	{
		float m_X;
		float m_Y;
		float m_Z;

		friend bool isEqual(const JBrainComponent& lhs, const JBrainComponent& rhs);		

		inline void setEqual(const JBrainComponent& rhs)
		{
			m_X = rhs.m_X;
			m_Y = rhs.m_Y;
			m_Z = rhs.m_Z;
		}

		inline void constrainLocation(const float& minX, const float& minY, const float& minZ,
			const float& maxX, const float& maxY, const float& maxZ)
		{
			m_X = fminf(maxX, fmaxf(minX, m_X));
			m_Y = fminf(maxY, fmaxf(minY, m_Y));
			m_Z = fminf(maxZ, fmaxf(minZ, m_Z));
		}

		JBrainComponent(float x, float y, float z) : m_X(x), m_Y(y), m_Z(z) {}
	};

	struct JDendrite : public JBrainComponent
	{
		// The weight it multiplies inputs by before passing the value
		// on to the body of the neuron:
		float m_weight;

		// The direction and distance to its most recent biggest neuron-fire contributor:
		float m_biggestInputX;
		float m_biggestInputY;
		float m_biggestInputZ;
		float m_biggestInputValue;
		float m_biggestInputDistance;
		bool m_biggestInputIsEnvironmentAxon;

		// The direction and distance to its nearest axon:
		float m_nearestAxonX;
		float m_nearestAxonY;
		float m_nearestAxonZ;
		bool m_nearestAxonPartOfSameNeuron;
		float m_nearestAxonDistance;
		bool m_nearestAxonIsEnvironmentAxon;
		
		// Its current output (multiplied by its weight).		
		float m_currentValue;
		bool m_recentInputAvailable;

		// As part of providing information to CGP updaters, we need to
		// track each dendrite's input to its neuron:
		std::vector<double> m_inputsProvided;
		double getAverageInput() const
		{
			// config.yaml variable name: INPUT_MAGNITUDE
			return std::accumulate(m_inputsProvided.begin(),
				m_inputsProvided.end(), 0) / static_cast<double>(m_inputsProvided.size());
		}

		void clearInputs() { m_inputsProvided.clear(); }

		JDendrite(float x, float y, float z, float weight):
			JBrainComponent(x, y, z), m_weight(weight),
			m_biggestInputX(0.0), m_biggestInputY(0.0), m_biggestInputZ(0.0),
			m_biggestInputValue(0.0), m_biggestInputDistance(0.0),
			m_biggestInputIsEnvironmentAxon(false), m_nearestAxonX(0.0),
			m_nearestAxonY(0.0), m_nearestAxonZ(0.0),
			m_nearestAxonPartOfSameNeuron(false),
			m_nearestAxonDistance(0.0), m_nearestAxonIsEnvironmentAxon(false),
			m_currentValue(-1.0), m_recentInputAvailable(false)
		{}

		// Equality operator:
		friend bool operator==(const JDendrite& lhs, const JDendrite& rhs);
		
	};

	struct JAxon : JBrainComponent
	{
		// The direction and distance to its nearest dendrite:
		float m_nearestDendriteX;
		float m_nearestDendriteY;
		float m_nearestDendriteZ;
		bool m_nearestDendriteIsPartOfSameNeuron;
		float m_nearestDendriteDistance;
		bool m_nearestDendriteIsActionDendrite;

		JAxon(float x, float y, float z) :
			JBrainComponent(x, y, z), m_nearestDendriteX(0.0),
			m_nearestDendriteY(0.0), m_nearestDendriteZ(0.0),
			m_nearestDendriteIsPartOfSameNeuron(false),
			m_nearestDendriteDistance(0.0),
			m_nearestDendriteIsActionDendrite(false)
		{}

		friend bool operator==(const JAxon& lhs, const JAxon& rhs);
	};

	struct JNeuron : JBrainComponent
	{
		std::vector<JAxon> m_axons;
		std::vector<JDendrite> m_dendrites;		
		
		// The output value the neuron creates when it fires:
		float m_fireValue;

		// The neuron's current health:
		float m_health;

		// Track how often we had the opportunity to fire because we weren't
		// on cooldown anymore and how often we did fire:
		unsigned int m_fireOpportunitiesSinceLastUpdate;
		unsigned int m_timesFiredSinceLastUpdate;

		// For tracking purposes:
		unsigned int m_neuronNumber;
		int m_timeStepsSinceLastFire;
		unsigned int m_age;

		friend bool operator==(const JNeuron& lhs, const JNeuron& rhs);

		JNeuron(const float& x, const float& y, const float& z,
			const float& fireValue, const float& health,
			const unsigned int& neuronNumber) :
			JBrainComponent(x, y, z), m_fireValue(fireValue),
			m_health(health), m_fireOpportunitiesSinceLastUpdate(0),
			m_timesFiredSinceLastUpdate(0), m_neuronNumber(neuronNumber),
			m_timeStepsSinceLastFire(-1), m_age(0)
		{}

		double getPercentageFire() const
		{
			double retVal = 0.0;
			// Make sure we don't divide by zero:
			if (m_fireOpportunitiesSinceLastUpdate != 0 &&
				m_timesFiredSinceLastUpdate != 0)
			{
				retVal = static_cast<double>(m_timesFiredSinceLastUpdate) /
					static_cast<double>(m_fireOpportunitiesSinceLastUpdate);
			}

			return retVal;			
		}

		void writeSelfHumanReadable(std::ostream& out) const
		{
			out << "### Begin neuron " << m_neuronNumber << " ###" << std::endl;
			out << "Body coordinates: (" << m_X << ", " << m_Y << ", " << m_Z << ")" << std::endl;
			out << "Axon coordinates: " << std::endl;
			for (unsigned int i = 0; i < m_axons.size(); ++i)
			{
				out << "\t" << i << ": (" << m_axons[i].m_X << ", "
					<< m_axons[i].m_Y << ", " << m_axons[i].m_Z
					<< ")" << std::endl;
			}

			out << "Dendrite coordinates: " << std::endl;
			for (unsigned int i = 0; i < m_dendrites.size(); ++i)
			{
				out << "\t" << i << ": (" << m_dendrites[i].m_X << ", "
					<< m_dendrites[i].m_Y << ", " << m_dendrites[i].m_Z
					<< ")" << std::endl;
			}
			out << "### End neuron " << m_neuronNumber << " ###" << std::endl;
		}
	};
}