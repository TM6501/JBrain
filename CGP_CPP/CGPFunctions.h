#pragma once

/*****************************************************************
* These are the functions that CGP neurons can choose from for
* their calculations. Any given CGP individual must choose the
* inputs and outputs properly as all must match.
* 
* These functions were taken from the original Python version of
* this library, CGPFunctions.py. Many unused functions were
* omitted.
*/

#include <map>
#include <functional>
#include <string>

namespace CGPFunctions
{
	class boolIn_boolOut
	{
	public:
		// P is often not available, give it a default:
		static bool And(bool X, bool Y, bool P = true);
		static bool Or(bool X, bool Y, bool P = true);
		static bool Nand(bool X, bool Y, bool P = true);
		static bool Nor(bool X, bool Y, bool P = true);
		static bool Xor(bool X, bool Y, bool P = true);
		static bool AndNotY(bool X, bool Y, bool P = true);

		// Return a pointer to a function based on a string describing it.
		// This allows the class to save off a list of strings and retriev
		// its same list of functions the next time that the program is run.
		static std::function<bool(bool, bool, bool)> getFuncFromString(const std::string& funcName);		
	};

	class doubleIn_doubleOut
	{
		// With double boolean functions, > 0 is considered "true" and
		// <= 0 is considered "false":
	public:
		static double And(double X, double Y, double P);  // X && Y
		static double Or(double X, double Y, double P);  // X || Y
		static double Nand(double X, double Y, double P);  // !(X && Y)
		static double Nor(double X, double Y, double P);  // (!(X || Y)
		static double Xor(double X, double Y, double P);  // (X || Y) && !(X && Y)
		static double AndNotY(double X, double Y, double P);  //X && !Y

		static double Add(double X, double Y, double P);  // X + Y
		static double Subtract(double X, double Y, double P);  // X - Y
		static double CMult(double X, double Y, double P); // X * P
		static double Mult(double X, double Y, double P);  // X * Y
		static double CDivide(double X, double Y, double P); // X / P
		static double Divide(double X, double Y, double P); // X / Y
		static double Inv(double X, double Y, double P);  // 1 / X. X if that throws an error.
		static double Abs(double X, double Y, double P);  // abs(X)
		static double SqrtX(double X, double Y, double P); // sqrt( abs(X) )
		static double SqrtXY(double X, double Y, double P);  // sqrt(X ^ 2 + Y ^ 2) / sqrt(2)
		static double CPow(double X, double Y, double P); // abs(X) ^ (P + 1)
		static double Pow(double X, double Y, double P);  // abs(X) ^ abs(Y)
		static double ExpX(double X, double Y, double P);  // ((e ^ X) - 1) / (e - 1)
		static double SinX(double X, double Y, double P);  // sin(X)
		static double ASinX(double X, double Y, double P);  // asin(X)
		static double CosX(double X, double Y, double P);  // cos(X)
		static double ACosX(double X, double Y, double P);  // acos(X)
		static double TanX(double X, double Y, double P);  // tan(X)
		static double ATanX(double X, double Y, double P);  // atan(X)
		static double LTE(double X, double Y, double P);  // double(X <= Y)  0 or 1
		static double GTE(double X, double Y, double P);  // double(X >= Y)  0 or 1
		static double GTEP(double X, double Y, double P);  // double(X >= P)  0 or 1
		static double LTEP(double X, double Y, double P);  // double(X <= P)  0 or 1
		static double Max(double X, double Y, double P);  // max(X, Y)
		static double Min(double X, double Y, double P);  // min(X, Y)
		static double XWire(double X, double Y, double P);  // X
		static double YWire(double X, double Y, double P);  // Y
		static double Const(double X, double Y, double P);  // P

		// Return a pointer to a function based on a string describing it.
		// This allows the class to save off a list of strings and retriev
		// its same list of functions the next time that the program is run.
		static std::function<double(double, double, double)> getFuncFromString(const std::string& funcName);
	};

	class neuronActivation
	{
		// These are functions designed to take in a single value
		// and return its value processed by a typical neural activation
		// function.  Typipcally, the input value would be a weighted
		// sum of inputs.
	public:
		static double Sigmoid(double X);
		static double Tanh(double X);
		static double Relu(double X);

		// Return a pointer to a function based on a string describing it.
		// This allows the class to save off a list of strings and retriev
		// its same list of functions the next time that the program is run.
		static std::function<double(double)> getFuncFromString(const std::string& funcName);
	};


}
