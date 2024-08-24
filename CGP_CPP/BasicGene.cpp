#include "pch.h"
#include "BasicGene.h"

extern unsigned int DEBUG_LEVEL;

namespace CGP
{
	BasicGene::BasicGene(GENE_TYPE type, int x, int y, int function, double p,
		double scaleVal)
		: type(type), X(x), Y(y), F(function), P(p), scale(scaleVal)
	{}
	BasicGene::~BasicGene() {}

	NoScaleGene::NoScaleGene(GENE_TYPE type, int x, int y, int function, double p)
		: type(type), X(x), Y(y), F(function), P(p)
	{}

	NoScaleGene::~NoScaleGene() {}
}
