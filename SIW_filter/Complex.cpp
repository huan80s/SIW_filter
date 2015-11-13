#include "Complex.h"
double complex::amplitude()
{
	double c;
	c = real*real + image*image;
	c = sqrt(c);
	return c;
}
complex complex::operator +(complex c2)//重载的实现
{
	complex c;
	c.real = c2.real + real;
	c.image = c2.image + image;
	return complex(c.real, c.image);
}
complex complex::operator -(complex c2)//重载的实现
{
	complex c;
	c.real = real - c2.real;
	c.image = image - c2.image;
	return complex(c.real, c.image);
}
complex complex::operator *(complex c2)
{
	complex c;
	c.real = real * c2.real - image * c2.image;
	c.image = real * c2.image + image * c2.real;
	return complex(c.real, c.image);
}

complex complex::operator *(double c2)
{
	complex c;
	c.real = real * c2;
	c.image = image * c2;
	return complex(c.real, c.image);
}

complex complex::operator /(complex c2)
{
	complex c;
	c.real = (real * c2.real + image * c2.image) / (c2.real * c2.real + c2.image * c2.image);
	c.image = (image * c2.real - real * c2.image) / (c2.real * c2.real + c2.image * c2.image);
	return complex(c.real, c.image);
}
void complex::display()
{
	cout << "(" << real << "," << image << ")" << endl;
}



