#ifndef Complex_h
#define Complex_h
#include <iostream>

using namespace std;

class complex//����������
{
private:
public:
	double real;
	double image;
	complex(double r = 0.0, double i = 0.0)//���캯��
	{
		real = r;
		image = i;
	}
	complex operator+(complex c2);//  +����Ϊ��Ա����
	complex operator-(complex c2);//  -����Ϊ��Ա����
	complex operator *(complex c2);//  *����Ϊ��Ա����
	complex operator* (double c2);
	complex operator /(complex c2);//  /����Ϊ��Ա����
	double amplitude();
	void display();

};

#endif // Complex_h
