#include <iostream>
#include "Complex.h"
#include <fstream>

void getB(complex* right, complex**A, int **Coef_Location, complex *x, int sizeM, int sizeN, int startpos)
{
	int N = 8;
	int n;
	complex zero;
	zero.real = 0; zero.image = 0;

	complex* temp = new complex[sizeM];
	for (int i = 0; i != sizeM; i++)
		temp[i] = zero;

	for (int i = 0; i != sizeM; i++)
	{
		for (int j = 0; j != N; j++)
		{
			if (Coef_Location[i][j] >= 0)
			{
				n = Coef_Location[i][j];
				temp[i] = temp[i] + A[i][j] * x[n];
			}
		}
	}

	for (int i = 0; i != sizeM; i++)
		right[i + startpos] = right[i + startpos] - temp[i];

	delete[] temp;
};