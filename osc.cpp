#include <iostream>

#include <cmath>

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <string>

#include "matrix.h"
#include "Approximation.h"

#include "osc.h"

using namespace std;

const double PI = 3.14159265359;

matrix count_freq(matrix A, int y_ind, int x_ind) /// dy/dx
{cout << "count_freq(matrix A, int y_ind, int x_ind)\n";A.info();

    matrix W(A.get_height(), 1);
    double h = A.get_data(1, x_ind) -  A.get_data(0, x_ind);

    W.set_data(0, 0, (-3*A.get_data(0, y_ind) + 4*A.get_data(1, y_ind) - A.get_data(2, y_ind))/2/h );
    for(int i = 1; i < A.get_height()-1; i++)
        W.set_data(i, 0, (A.get_data(i+1, y_ind) - A.get_data(i-1, y_ind))/2/h );
    W.set_data(W.get_height()-1, 0, (A.get_data(A.get_height()-3, y_ind) - 4*A.get_data(W.get_height()-2, y_ind) + 3*A.get_data(W.get_height()-1, y_ind))/2/h );

    return W;
};

matrix count_energy(matrix A, double mz_a)
{
    double Ek, Em, Ea, a = PI/180.0;
    matrix E(A.get_height(), 3);

    for(int i = 0; i < A.get_height(); i++)
    {
        Ek = A.get_data(i, 2)*a*A.get_data(i, 2)*a/2.0*0.0001; ///*I
        Em = A.get_data(i, 1)*a*mz_a;
        Ea = Ek+Em;

        E.set_data(i, 0, Ek);
        E.set_data(i, 1, Em);
        E.set_data(i, 2, Ea);
    }
    return E;
}

