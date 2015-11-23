#include "Complex.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
void show_binc(int k, complex* binc, int num_total){
	char filename[256];
	sprintf_s(filename,"binc\\binc%d.txt", k);
	ofstream out(filename);
	for (int i = 0; i != num_total; i++){
		if (binc[i].real != 0 || binc[i].image != 0)
			out << binc[i].real << " " << binc[i].image << " " << i << endl;
	}
	out.close();
}