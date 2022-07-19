#ifndef COMPLEX_H
#define COMPLEX_H

template<typename T>
class Complex // Complex number
{
public:
    T r, i;
    Complex();
    Complex(T r, T i);
    Complex operator+(const Complex& c) const;
    Complex operator-(const Complex& c) const;
    void operator+=(const Complex& c);
    void operator-=(const Complex& c);
    void operator*=(const Complex& c);
    Complex operator/(float s);


    void display();
    T magnitude_squared();
    T magnitude();
};

#endif // COMPLEX_H
