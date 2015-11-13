#ifndef BASIC_FUNC_H
#define BASIC_FUNC_H
#include "Complex.h"

void inner_5node(int j, int i, int i_1, int i_2, int i_3, int i_4, complex** Coef, int** Coef_location);
void inner_4node(int j, int i, int i_1, int i_2, int i_3, complex** Coef, int** Coef_location);
void boundary_4node(int j, int i, int i_1, int i_2, int i_3, complex** Coef, int** Coef_location);
void corner_4node(int j, int i, int i_1, int i_2, int i_3, complex** Coef, int** Coef_location);

#endif