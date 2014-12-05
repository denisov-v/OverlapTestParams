#include "biglong.h"

//stores all computed values for the reccurence formula
class TValueMatrix
{
	blong::biglong** T_;
	int K_;
	int M_;
	
	TValueMatrix() = delete;
	TValueMatrix(const TValueMatrix&) = delete;
	TValueMatrix(const TValueMatrix&&) = delete;
	TValueMatrix& operator=(const TValueMatrix&) = delete;
	TValueMatrix& operator=(const TValueMatrix&&) = delete;
public:

	TValueMatrix(int K, int M) : K_(K), M_(M)
	{
		T_ = new blong::biglong*[K_];
		for (int i = 0; i < K_; ++i)
			T_[i] = new blong::biglong[M_ + 2];
	}	

	~TValueMatrix()
	{
		for (int i = 0; i < K_; ++i)
			delete[] T_[i];
		delete[] T_;
	}	

	blong::biglong& operator()(int alpha, int i)
	{
		return T_[alpha][i + 1];
	}
};

//implements the reccurence formula for occurrence probabilities of the template
//in The Overlapping Template Matching Test from NIST Statistical Test Suite
//parameters
//*IN* K: number of degrees of freedom
//*IN* M: length in bits of a string to be tested
//*IN* m: length in bits of the template
//*OUT* pi: occurrence probabilities of the template
inline void compute_probabilities(int K, int M, int m, double* pi)
{
	TValueMatrix T(K, M);

	T(0, -1) = 1;
	T(0, 0) = 1;
	for (int i = 1; i <= m - 1; ++i)
		T(0, i) = blong::biglong::two * T(0, i - 1);
	for (int i = m; i <= M; ++i)
		T(0, i) = (blong::biglong::two * T(0, i - 1)).trunc_sub(T(0, i - m - 1));

	for (int i = -1; i <= m - 1; ++i)
		T(1, i) = 0;
	T(1, m) = 1;
	T(1, m + 1) = 2;
	for (int i = m + 2; i <= M; ++i)
	{
		T(1, i) = 0;
		for (int j = -1; j <= i - m - 1; ++j)
			T(1, i) += T(0, j) * T(0, i - m - 2 - j);
	}

	for (int alpha = 2; alpha < K; ++alpha)
	{
		T(alpha, -1) = 0;
		for (int i = 0; i <= M; ++i)
		{
			T(alpha, i) = T(alpha - 1, i - 1);
			for (int j = -1; j <= i - 2 * m - alpha; ++j)
				T(alpha, i) += T(0, j) * T(alpha - 1, i - m - 2 - j);
		}
	}

	pi[0] = T(0, M).bit_rshift_frac_part(M);
	pi[1] = T(1, M).bit_rshift_frac_part(M);
	pi[2] = T(2, M).bit_rshift_frac_part(M);
	pi[3] = T(3, M).bit_rshift_frac_part(M);
	pi[4] = T(4, M).bit_rshift_frac_part(M);
	pi[5] = 1 - pi[0] - pi[1] - pi[2] - pi[3] - pi[4];
}