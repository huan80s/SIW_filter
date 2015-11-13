#include "Complex.h"
#include <iostream>
#include <fstream>


void toEvalue(complex *x, int num_total)
{
	ofstream out_e("Evalue.txt");
	for (int i = 0; i != num_total; i++)
		out_e << x[i].amplitude() << endl;
	out_e.close();
}
