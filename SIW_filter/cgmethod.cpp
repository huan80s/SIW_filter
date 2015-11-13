#include <iostream>
#include <fstream>
#include "Complex.h"
#include "Vector.h"
#include <vector>
using std::vector;
using std::cout;

//投影分解法使用的测试方法
void cgmethod(complex** A, int** Coef_Location, complex* right, complex* x, int sizeM, int sizeN, int N_matrix)
{
	//cout << "test cgmethod!" << endl;

	complex zero;
	zero.real = 0;
	zero.image = 0;

	int n;	//位置变量
	//*******************
	complex* r = new complex[sizeN];
	complex* p = new complex[sizeN];
	complex* p_conj = new complex[sizeN];
	complex* r_conj = new complex[sizeN];
	complex* Ap = new complex[sizeN];
	complex* Ap_conj = new complex[sizeN];
	complex* B = new complex[sizeN];

	complex* b = new complex[sizeN];
	for (int i = 0; i != sizeN; i++)
		b[i] = zero;
	//b=AH*right
	/*ofstream out_right("test_rightBinc.txt");
	for (int i = 0; i != sizeM; i++)
	out_right << right[i] << endl;
	out_right.close();*/
	for (int i = 0; i != sizeM; i++)
	{
		for (int j = 0; j != N_matrix; j++)
		{
			if (Coef_Location[i][j] >= 0)
			{
				n = Coef_Location[i][j];		//	位置修正
				b[n] = b[n] + conj(A[i][j])*right[i];
			}
		}
	}
	/*ofstream out_b("test_Ahb.txt");
	for (int i = 0; i != sizeN; i++)
	out_b << b[i] << endl;
	out_b.close();*/

	//****************************************
	for (int i = 0; i != sizeN; i++)
	{
		x[i].real = 0.0;
		x[i].image = 0.0;
	}
	for (int i = 0; i != sizeN; i++)
	{
		B[i] = b[i];
		r[i] = b[i];
		p[i] = b[i];
		p_conj[i] = conj(b[i]);
		r_conj[i] = conj(b[i]);
	}
	//double cond = VectorNorm(sizeN, r) / VectorNorm(sizeN, B);
	double cond = 1;
	double fenzi = 0.0;
	double fenmu = 0.0;
	complex a;
	double eps = 1e-10;
	int time = 0;
	complex* Ap_temp = new complex[sizeM];		//储存AH*(A*P)的中间变量
	complex* Ap_conj_temp = new complex[sizeM];

	//*****************************************************************
	while (cond>eps)
	{
		//cout << "cond...." << cond << endl;
		//Ap 要求Ap 实际上求的是AH*A*p; (N*M)(M*N)*(N*1)=(N*M)(M*1)=(N*1)
		//A(M*N)*p(N*1)
		for (int i = 0; i != sizeN; i++)
			Ap[i] = zero;
		for (int i = 0; i != sizeM; i++)
			Ap_temp[i] = zero;

		for (int i = 0; i != sizeM; i++)
		{
			for (int j = 0; j != N_matrix; j++)
			{
				if (Coef_Location[i][j] >= 0)
				{
					n = Coef_Location[i][j];         //修正位置
					Ap_temp[i] = Ap_temp[i] + A[i][j] * p[n];
				}
			}
		}
		//AH*Ap_temp
		for (int i = 0; i != sizeM; i++)
		{
			for (int j = 0; j != N_matrix; j++)
			{
				if (Coef_Location[i][j] >= 0)
				{
					n = Coef_Location[i][j];
					Ap[n] = Ap[n] + conj(A[i][j])*Ap_temp[i];
				}
			}
		}
		/*ofstream out_Ap("test_Ap.txt");
		for (int i = 0; i != sizeN; i++)
		out_Ap << Ap[i] << endl;
		out_Ap.close();*/

		//**************************************A*p
		complex ak = Vectorneiji(sizeN, r_conj, r) / Vectorneiji(sizeN, p_conj, Ap);
		/*cout << Vectorneiji(sizeN, r_conj, r)<<endl;
		cout << Vectorneiji(sizeN, p_conj, Ap) << endl;
		cout << Vectorneiji(sizeN, r_conj, r) / Vectorneiji(sizeN, p_conj, Ap) << endl;
		cout << "ak:" << ak << endl;*/
		fenzi = 0.0;
		fenmu = 0.0;
		for (int i = 0; i != sizeN; i++)
		{
			a = ak*p[i];
			x[i] = x[i] + a;
			r[i] = r[i] - ak*Ap[i];
			fenzi = fenzi + a.real*a.real + a.image*a.image;
			fenmu = fenmu + x[i].real*x[i].real + x[i].image*x[i].image;
		}
		cond = sqrt(fenzi / fenmu);

		//Ap_conj
		//A*p_conj
		for (int i = 0; i != sizeN; i++)
			Ap_conj[i] = zero;
		for (int i = 0; i != sizeM; i++)
			Ap_conj_temp[i] = zero;

		for (int i = 0; i != sizeM; i++)
		{
			for (int j = 0; j != N_matrix; j++)
			{
				if (Coef_Location[i][j] >= 0)
				{
					n = Coef_Location[i][j];         //修正位置
					Ap_conj_temp[i] = Ap_conj_temp[i] + A[i][j] * p_conj[n];
				}
			}
		}
		//AH*Ap_temp
		for (int i = 0; i != sizeM; i++)
		{
			for (int j = 0; j != N_matrix; j++)
			{
				if (Coef_Location[i][j] >= 0)
				{
					n = Coef_Location[i][j];
					Ap_conj[n] = Ap_conj[n] + conj(A[i][j])*Ap_conj_temp[i];
				}
			}
		}
		/*ofstream out_Ap_conj("test_Ap_conj.txt");
		for (int i = 0; i != sizeN; i++)
		out_Ap_conj << Ap_conj[i] << endl;
		out_Ap_conj.close();*/

		//******************************
		for (int i = 0; i != sizeN; i++)
		{
			r_conj[i] = r_conj[i] - conj(ak)*Ap_conj[i];
		}
		complex minus1(-1, 0);
		complex bk = minus1*Vectorneiji(sizeN, Ap_conj, r) / Vectorneiji(sizeN, p_conj, Ap);
		for (int kk = 0; kk != sizeN; kk++)
		{
			p[kk] = r[kk] + bk*p[kk];
			p_conj[kk] = r_conj[kk] + conj(bk)*p_conj[kk];
		}
		time++;
	}//while
	/*cout << "迭代次数：" << time << endl;
	cout << "end cg_method." << endl;*/


	delete[] r;
	delete[] p;
	delete[] p_conj;
	delete[] r_conj;
	delete[] Ap;
	delete[] Ap_conj;
	delete[] B;
};

////共轭梯度子函数 求解A*x=b
////说明：A为系数矩阵
////      b为右端向量
////      size为矩阵阶数
////       输出为解 x
////*******************************************************************
//complex* cg_method(complex** A, complex* q, int size, int disp, int** Coef_location)
//{
//	complex* x1 = new complex[size];
//	complex* r = new complex[size];
//	complex* p = new complex[size];
//	complex* p_conj = new complex[size];
//	complex* r_conj = new complex[size];
//	complex* Ap = new complex[size];
//	complex* Ap_conj = new complex[size];
//	complex* B = new complex[size];
//	//****************************************
//	int k = 0, i = 0, j = 0, n = 0;
//
//	for (k = 0; k != size; k++)
//	{
//		x1[k].real = 0.0;
//		x1[k].image = 0.0;
//	}
//	for (k = 0; k != size; k++)
//	{
//		B[k] = q[k];
//		r[k] = q[k];
//		p[k] = q[k];
//		p_conj[k] = conj(q[k]);
//		r_conj[k] = conj(q[k]);
//	}
//	double cond = 1;
//	double eps = 1e-10;
//	//*****************************************************************
//	int mm = 0;
//	while (cond>eps)
//	{
//		for (k = 0; k != size; k++)
//		{
//			Ap[k].real = 0.0;
//			Ap[k].image = 0.0;
//		}
//		for (k = 0; k != size; k++)
//		{
//			for (j = 0; j != N_matrix; j++)
//			{
//				if (Coef_location[k][j] >= 0)
//				{
//					n = Coef_location[k][j];
//					Ap[k] = Ap[k] + A[k][j] * p[n];
//				}
//			}
//			//   Ap[k] = Ap[k] + A[k][j]*p[j];
//		}
//		//**************************************A*p
//		complex ak = Vectorneiji(size, r_conj, r) / Vectorneiji(size, p_conj, Ap);
//		for (i = 0; i != size; i++)
//		{
//			x1[i] = x1[i] + ak*p[i];;
//			r[i] = r[i] - ak*Ap[i];
//		}
//		for (k = 0; k != size; k++)
//		{
//			Ap_conj[k].real = 0.0;
//			Ap_conj[k].image = 0.0;
//		}
//		for (k = 0; k != size; k++)
//		{
//			for (j = 0; j != N_matrix; j++)
//			{
//				if (Coef_location[k][j] >= 0)
//				{
//					n = Coef_location[k][j];
//					Ap_conj[n] = Ap_conj[n] + conj(A[k][j]) * p_conj[k];
//				}
//			}//Ap_conj[k] = Ap_conj[k] + conj(A[j][k])*p_conj[j];
//		}
//		for (k = 0; k != size; k++)
//		{
//			r_conj[k] = r_conj[k] - conj(ak)*Ap_conj[k];
//		}
//		complex temp11(-1, 0);
//		complex bk = temp11*Vectorneiji(size, Ap_conj, r) / Vectorneiji(size, p_conj, Ap);
//		for (k = 0; k != size; k++)
//		{
//			p[k] = r[k] + bk*p[k];
//			p_conj[k] = r_conj[k] + conj(bk)*p_conj[k];
//		}
//		cond = VectorNorm(size, r) / VectorNorm(size, B);
//		if (disp == 1)//disp=1时，显示每次迭代的比较结果
//			cout << cond << endl;
//		mm++;
//	}//while
//	cout << "迭代次数：" << mm << endl;
//	delete[] r;
//	delete[] p;
//	delete[] p_conj;
//	delete[] r_conj;
//	delete[] Ap;
//	delete[] Ap_conj;
//	delete[] B;
//	return x1;
//}//cg_method
