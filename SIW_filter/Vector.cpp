/****************** �����һЩ������C=A*B;b=A*b;****************************/
/****************** ���������С���˷�At*A*x=At*b; ****************************/
#include <cmath>
#include "Complex.h"
#include "Vector.h"
#include <vector>

using std::vector;
//����************************************************************
complex conj(complex A)
{
	complex t(0, 0);
	t.real = A.real;
	t.image = -A.image;
	return t;
};
double VectorNorm(int const num,complex*const B)                            // ||a||_2;ģֵ
{
	double a = 0;
	for(int i=0;i<num;i++)
		a = a + B[i].real*B[i].real + B[i].image*B[i].image;
	a = sqrt(a);
	return a;
}

double VectorNorm(double*const b,int const num)                             //ʵ����ģֵ
{
    double temp=0;
    for(int i=0;i!=num;i++)
    {
    temp+=b[i]*b[i];
    }
    temp=sqrt(temp);
    return temp;
}

double VectorNorm(complex*const B,int const num)
{
	double a = 0;
	for(int i=0;i<num;i++)
		a = a + B[i].real*B[i].real + B[i].image*B[i].image;
	a = sqrt(a);
	return a;
}


complex Vectorneiji(int const num,complex*const A,complex*const B)           //����a[]*b[];
{
	complex t;
	t.real=0;
	t.image=0;
	for(int i=0;i<num;i++)
	{
		t.real=t.real+A[i].real*B[i].real+A[i].image*B[i].image;
		t.image=t.image +A[i].real*B[i].image-A[i].image*B[i].real;
	}
	return  t;

}

void ConjMatrix(complex**const a,complex** b,int const m,int const n)        //Ah;����ת��; b(M*N)
{
    for(int i=0;i!=n;i++)
        for(int j=0;j!=m;j++)
        {
            b[j][i]=conj(a[i][j]);
        }
}

void ConjMatrix(complex**const a,complex** b,int const m)                    //Ah;����Ĺ���ת��
{
    for(int i=0;i!=m;i++)
        for(int j=0;j!=m;j++)
            b[j][i]=conj(a[i][j]);
}

complex* ConjVector(complex*const a,int const m)                            //bh;�����Ĺ���;
{
    complex* b;
    b=new complex[m];
    for(int j=0;j!=m;j++)
            b[j]=conj(a[j]);
    return b;
}

void Inner_product(complex** a,complex* b,complex* c,int m,int n)           //����c=A*b;
{
	for(int i=0;i!=m;i++)
    {
        for(int j=0;j!=n;j++)
            c[i]=c[i]+a[i][j]*b[j];
    }
}

void Inner_product(complex** a, vector<complex>& b, complex* c, int m, int n)           //����c=A*b;
{
	for (int i = 0; i != m; i++)
	{
		for (int j = 0; j != n; j++)
			c[i] = c[i] + a[i][j] * b[j];
	}
}


double Inner_product(double *a,double* b,int n)                             //ʵ�������˻�
{
    double temp=0;
    for(int i=0;i!=n;i++)
    {
        temp=temp+a[i]*b[i];
    }
    return temp;
}

double Inner_product(complex* B,int n)
{
    double a = 0;
	for (int i = 0; i < n; i++)
	{
		a =a+ B[i].real*B[i].real + B[i].image*B[i].image;
	}
	return a;
}

void MulMatrix(complex**const a,complex**const b,complex **C,int const m,int const n)     //����C=A*B
{
   int i,j,k;
   for(i=0;i<m;i++)
   {
       for(j=0;j<m;j++){
         for(k=0;k<n;k++)
            C[i][j]=C[i][j]+a[i][k]*b[k][j];
     }
   }

}

void MulMatrix(complex**const a,complex**const b,complex **C,int const m)     //����C=A*B
{
   int i,j,k;
   for(i=0;i<m;i++)
   {
       for(j=0;j<m;j++){
         for(k=0;k<m;k++)
            C[i][j]=C[i][j]+a[i][k]*b[k][j];
     }
   }

}
