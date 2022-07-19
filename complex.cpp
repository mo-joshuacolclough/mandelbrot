#include <iostream>
#include "complex.h"
#include <cmath>

template<typename T>
Complex<T>::Complex(){}

template<typename T>
Complex<T>::Complex(T r_, T i_)
{
    r = r_;
    i = i_;
}

template<typename T>
void Complex<T>::display()
{
    std::cout << r << " + " << i << "i" << std::endl;
}

template<typename T>
Complex<T> Complex<T>::operator+(const Complex<T>& c) const
{
    Complex summed = Complex(0, 0);
    summed.r = r + c.r;
    summed.i = i + c.i;

    return summed;
}

template<typename T>
Complex<T> Complex<T>::operator-(const Complex<T>& c) const
{
    return Complex(r - c.r, i - c.i);
}

template<typename T>
Complex<T> Complex<T>::operator/(float s)
{
    return Complex<T>(r/s, i/s);
}

template<typename T>
void Complex<T>::operator+=(const Complex<T>& c)
{
    r += c.r;
    i += c.i;
}

template<typename T>
void Complex<T>::operator-=(const Complex<T>& c)
{
    r -= c.r;
    i -= c.i;
}

template<typename T>
void Complex<T>::operator*=(const Complex<T>& c)
{
    r *= c.r;
    i *= c.i;
}


template<typename T>
T Complex<T>::magnitude_squared()
{
    return r * r + i * i;
}

template<typename T>
T Complex<T>::magnitude()
{
    return std::sqrt(magnitude_squared());
}

//template class Complex<long double>;
//template class Complex<double>;
//template class Complex<float>;
//template class Complex<int>;

