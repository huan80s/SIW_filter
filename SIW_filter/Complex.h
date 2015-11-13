#ifndef Complex_h
#define Complex_h
#include <iostream>

using namespace std;

class complex//复数类声明
{
private:
public:
	double real;
	double image;
	complex(double r = 0.0, double i = 0.0)//构造函数
	{
		real = r;
		image = i;
	}
	complex operator+(complex c2);//  +重载为成员函数
	complex operator-(complex c2);//  -重载为成员函数
	complex operator *(complex c2);//  *重载为成员函数
	complex operator* (double c2);
	complex operator /(complex c2);//  /重载为成员函数
	double amplitude();
	void display();

};

#endif // Complex_h
