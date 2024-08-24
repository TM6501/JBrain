#pragma once
#include "Enums.h"

namespace CGP
{
	struct BasicGene
	{
		GENE_TYPE type;
		int X;  // Node number for first input. For output genes, this also serves
		        // as the selector for the node from which to pull the output.
		int Y;  // Node number for second input.
		int F;  // Index into list of functions.
		double P;  // P (weight in some cases) to apply to function.
		double scale;  // Scale to apply to final result.

		BasicGene(GENE_TYPE type, int x, int y, int function, double p,
			double scaleVal);
		~BasicGene();
	};

	struct NoScaleGene
	{
		GENE_TYPE type;
		int X;	
		int Y;  
		int F;
		double P;		

		NoScaleGene(GENE_TYPE type, int x, int y, int function, double p);
		~NoScaleGene();

		// Copy operator:
		NoScaleGene& operator=(const NoScaleGene& other)
		{
			if (this == &other)
				return *this;

			type = other.type;
			X = other.X;
			Y = other.Y;
			F = other.F;
			P = other.P;

			return *this;
		}

		friend bool operator==(const NoScaleGene& lhs, const NoScaleGene& rhs);
	};	
}
