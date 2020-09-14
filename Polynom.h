#ifndef POLYNOM_H_INCLUDED
#define POLYNOM_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

class polynom
{
private:

protected:
    int degree;
    double* coefficient;

public:
    ///constructors
    polynom();
    polynom(int n);
    polynom(double*, int n);
    ~polynom();
    polynom(const polynom&);


    ///
    int get_degree();
    double get_coeff(int index_degree);
    double* get_arr_coeff(int&);/// и степень в аргументе
    void set_degree();
    void set_coeff(double coeff, int index_degree);
    void set_coeff(double* , int n);


    void info();

};

#endif // POLYNOM_H_INCLUDED
