#include "pch.h"
#include "framework.h"
#include "CGPFunctions.h"
#include <cmath>
#include <string>
#include <algorithm>

namespace CGPFunctions
{
	bool boolIn_boolOut::And(bool X, bool Y, bool P)
	{
		return X && Y;
	}

	bool boolIn_boolOut::Or(bool X, bool Y, bool P)
	{
		return X || Y;
	}

	bool boolIn_boolOut::Nand(bool X, bool Y, bool P)
	{
		return !(X && Y);
	}

	bool boolIn_boolOut::Nor(bool X, bool Y, bool P)
	{
		return !(X || Y);
	}

	bool boolIn_boolOut::Xor(bool X, bool Y, bool P)
	{
		return ((X || Y) && !(X && Y));
	}

	bool boolIn_boolOut::AndNotY(bool X, bool Y, bool)
	{
		return X && !Y;
	}

	std::function<bool(bool, bool, bool)> boolIn_boolOut::getFuncFromString(const std::string& funcName)
	{
		std::string func = funcName;
		std::transform(func.begin(), func.end(), func.begin(), ::toupper);

		if (func == "AND")
			return boolIn_boolOut::And;
		else if (func == "OR")
			return boolIn_boolOut::Or;
		else if (func == "NAND")
			return boolIn_boolOut::Nand;
		else if (func == "NOR")
			return boolIn_boolOut::Nor;
		else if (func == "XOR")
			return boolIn_boolOut::Xor;
		else if (func == "ANDNOTY")
			return boolIn_boolOut::AndNotY;
		else
			return nullptr;
	}

	double doubleIn_doubleOut::And(double X, double Y, double P)
	{
		double retValue = -1.0;
		if (X > 0.0 && Y > 0.0)
			retValue = 1.0;
		
		return retValue;
	}

	double doubleIn_doubleOut::Or(double X, double Y, double P)
	{
		double retValue = -1.0;
		if (X > 0.0 || Y > 0.0)
			retValue = 1.0;

		return retValue;
	}

	double doubleIn_doubleOut::Nand(double X, double Y, double P)
	{
		return -And(X, Y, P);
	}

	double doubleIn_doubleOut::Nor(double X, double Y, double P)
	{
		return -Or(X, Y, P);
	}

	double doubleIn_doubleOut::Xor(double X, double Y, double P)
	{
		double retValue = -1.0;
		if ((X > 0.0 && Y <= 0.0) || (X <= 0.0 && Y > 0.0))
			retValue = 1.0;

		return retValue;
	}

	double doubleIn_doubleOut::AndNotY(double X, double Y, double P)
	{
		double retValue = -1.0;
		if (X > 0.0 and Y <= 0.0)
			retValue = 1.0;

		return retValue;
	}

	double doubleIn_doubleOut::Add(double X, double Y, double P)
	{
		return X + Y;
	}
	
	double doubleIn_doubleOut::Subtract(double X, double Y, double P)
	{
		return X - Y;
	}
	
	double doubleIn_doubleOut::CMult(double X, double Y, double P)
	{
		return X * P;
	}
	
	double doubleIn_doubleOut::Mult(double X, double Y, double P)
	{
		return X * Y;
	}

	double doubleIn_doubleOut::CDivide(double X, double Y, double P)
	{
		// Only divide if the divisor isn't close to 0:
		if (abs(P) < 0.001)
			return X;
		else
			return X / P;
	}
	
	double doubleIn_doubleOut::Divide(double X, double Y, double P)
	{
		// Only divide if the divisor isn't close to 0:
		if (abs(Y) < 0.001)
			return X;
		else
			return X / Y;
	}
	
	double doubleIn_doubleOut::Inv(double X, double Y, double P)
	{
		// Prevent division by 0 errors:
		double retValue = X;
		if (std::abs(X) > 0.0001)
			retValue = 1 / X;

		return retValue;
	}
	
	double doubleIn_doubleOut::Abs(double X, double Y, double P)
	{
		return std::abs(X);
	}

	double doubleIn_doubleOut::SqrtX(double X, double Y, double P)
	{
		return std::sqrt(std::abs(X));
	}

	double doubleIn_doubleOut::SqrtXY(double X, double Y, double P)
	{
		return double(std::sqrt((X * X) + (Y * Y)) / std::sqrt(2.0));
	}

	double doubleIn_doubleOut::CPow(double X, double Y, double P)
	{
		return std::pow(X, P);
	}
	
	double doubleIn_doubleOut::Pow(double X, double Y, double P)
	{
		return std::pow(std::abs(X), std::abs(Y));
	}

	double doubleIn_doubleOut::ExpX(double X, double Y, double P)
	{
		return double((std::exp(X) - 1.0) / (std::exp(1.0) - 1));
	}
	
	double doubleIn_doubleOut::SinX(double X, double Y, double P)
	{
		return std::sin(X);
	}

	double doubleIn_doubleOut::ASinX(double X, double Y, double P)
	{
		return std::asin(X);
	}
	
	double doubleIn_doubleOut::CosX(double X, double Y, double P)
	{
		return std::cos(X);
	}
	
	double doubleIn_doubleOut::ACosX(double X, double Y, double P)
	{
		return std::acos(X);
	}
	
	double doubleIn_doubleOut::TanX(double X, double Y, double P)
	{
		return std::tan(X);	
	}
	
	double doubleIn_doubleOut::ATanX(double X, double Y, double P)
	{
		return std::atan(X);	
	}
	
	double doubleIn_doubleOut::LTE(double X, double Y, double P)
	{
		return double(X <= Y);
	}
	
	double doubleIn_doubleOut::GTE(double X, double Y, double P)
	{
		return double(X >= Y);
	}

	double doubleIn_doubleOut::GTEP(double X, double Y, double P)
	{
		return double(X >= P);
	}
	
	double doubleIn_doubleOut::LTEP(double X, double Y, double P)
	{
		return double(X <= P);
	
	}
	
	double doubleIn_doubleOut::Max(double X, double Y, double P)
	{
		return std::fmax(X, Y);
	}

	double doubleIn_doubleOut::Min(double X, double Y, double P)
	{
		return std::fmin(X, Y);
	}
	
	double doubleIn_doubleOut::XWire(double X, double Y, double P)
	{
		return X;
	}

	double doubleIn_doubleOut::YWire(double X, double Y, double P)
	{
		return Y;
	}

	double doubleIn_doubleOut::Const(double X, double Y, double P)
	{
		return P;
	}

	std::function<double(double, double, double)> doubleIn_doubleOut::getFuncFromString(const std::string& funcName)
	{
		std::string func = funcName;
		std::transform(func.begin(), func.end(), func.begin(), ::toupper);

		if (func == "AND")
			return doubleIn_doubleOut::And;		
		else if (func == "OR")
			return doubleIn_doubleOut::Or;
		else if (func == "NAND")
			return doubleIn_doubleOut::Nand;
		else if (func == "NOR")
			return doubleIn_doubleOut::Nor;
		else if (func == "XOR")
			return doubleIn_doubleOut::Xor;
		else if (func == "ANDNOTY")
			return doubleIn_doubleOut::AndNotY;
		else if (func == "ADD")
			return doubleIn_doubleOut::Add;
		else if (func == "SUBTRACT")
			return doubleIn_doubleOut::Subtract;
		else if (func == "CMULT")
			return doubleIn_doubleOut::CMult;
		else if (func == "MULT")
			return doubleIn_doubleOut::Mult;
		else if (func == "CDIVIDE")
			return doubleIn_doubleOut::CDivide;
		else if (func == "DIVIDE")
			return doubleIn_doubleOut::Divide;
		else if (func == "INV")
			return doubleIn_doubleOut::Inv;
		else if (func == "ABS")
			return doubleIn_doubleOut::Abs;
		else if (func == "SQRTX")
			return doubleIn_doubleOut::SqrtX;
		else if (func == "SQRTXY")
			return doubleIn_doubleOut::SqrtXY;
		else if (func == "CPOW")
			return doubleIn_doubleOut::CPow;
		else if (func == "POW")
			return doubleIn_doubleOut::Pow;
		else if (func == "EXPX")
			return doubleIn_doubleOut::ExpX;
		else if (func == "SINX")
			return doubleIn_doubleOut::SinX;
		else if (func == "ASINX")
			return doubleIn_doubleOut::ASinX;
		else if (func == "COSX")
			return doubleIn_doubleOut::CosX;
		else if (func == "ACOSX")
			return doubleIn_doubleOut::ACosX;		
		else if (func == "TANX")
			return doubleIn_doubleOut::TanX;
		else if (func == "ATANX")
			return doubleIn_doubleOut::ATanX;
		else if (func == "LTE")
			return doubleIn_doubleOut::LTE;
		else if (func == "GTE")
			return doubleIn_doubleOut::GTE;
		else if (func == "GTEP")
			return doubleIn_doubleOut::GTEP;
		else if (func == "LTEP")
			return doubleIn_doubleOut::LTEP;
		else if (func == "MAX")
			return doubleIn_doubleOut::Max;
		else if (func == "MIN")
			return doubleIn_doubleOut::Min;
		else if (func == "XWIRE")
			return doubleIn_doubleOut::XWire;
		else if (func == "YWIRE")
			return doubleIn_doubleOut::YWire;
		else if (func == "CONST")
			return doubleIn_doubleOut::Const;
		else
			return nullptr;
	}

	double neuronActivation::Sigmoid(double X)
	{
		return double(1.0 / (1.0 + std::exp(-X)));
	}

	double neuronActivation::Tanh(double X)
	{
		return tanh(X);
	}

	double neuronActivation::Relu(double X)
	{
		double retValue = 0.0;
		if (X > 0.0)
			retValue = X;
		
		return retValue;
	}

	std::function<double(double)> neuronActivation::getFuncFromString(const std::string& funcName)
	{
		std::string func = funcName;
		std::transform(func.begin(), func.end(), func.begin(), ::toupper);

		if (func == "SIGMOID")
			return neuronActivation::Sigmoid;
		else if (func == "TANH")
			return neuronActivation::Tanh;
		else if (func == "RELU")
			return neuronActivation::Relu;
		else
			return nullptr;
	}
}