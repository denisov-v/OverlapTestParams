//this program computes occurrence probabilities of the template
//in The Overlapping Template Matching Test from NIST Statistical Test Suite

#include "prob.h"
#include <iomanip>

#define K 5 //number of degrees of freedom
#define M 1032 //length in bits of a string to be tested
#define m 9 //length in bits of the template
#define PRECISION 10 //decimal precision to format floating-point values

void main()
{
	//occurrence probabilities of the template
	double pi[K+1];

	compute_probabilities(K, M, m, pi);

	std::cout << std::setprecision(PRECISION);
	std::cout << pi[0] << std::endl;
	std::cout << pi[1] << std::endl;
	std::cout << pi[2] << std::endl;
	std::cout << pi[3] << std::endl;
	std::cout << pi[4] << std::endl;
	std::cout << pi[5] << std::endl;
}