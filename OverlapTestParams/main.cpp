#include "prob.h"
#include <iomanip>

#define K 5
#define M 1032
#define m 9
#define PRECISION 8

void main()
{
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