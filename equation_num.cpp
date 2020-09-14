#include <iostream>
#include <stdlib.h> ///rand()
#include "time.h"
#include <string>
#include <functional>


#include "cmath"
#include "Matrix.h"
#include "Approximation.h"
#include <fstream>

#include "equation_num.h"
using namespace std;
const double Pi = 3.1415926535897932384626433832795;
const double PI = 3.1415926535897932384626433832795;

/*DifEq_2order::DifEq_2order()
{
    fa0 = fa2 = fa4 = ga1 = ga3 = ga5 = 0.0;
    matrix nm;
        result = nm;
}*/

void base_equation(string Data_Way, string file_name, double b, double c)
{
    ofstream eq_fout(Data_Way + file_name + ".txt");
        cout << "equation: d2y/dx2  + " << b << "*dy/dx + "<< c << "*x = 0;\n";
        eq_fout << "equation: d2y/dx2  + " << b << "*dy/dx + "<< c << "*x = 0;\n";
    double D2 = 0;
    D2 = b*b - 4*c;

    if(D2 < 0)
    {
        double w = sqrt(D2*(-1.0))/2.0;
            cout << "exp(" << -1*b/2 << ")[Acos(" << w << "t) + Bsin(" << w << "t)]" << endl << "f = " << w/2.0/Pi << " Hz"<< endl;
         eq_fout << "exp(" << -1*b/2 << ")[Acos(" << w << "t) + Bsin(" << w << "t)]" << endl << "f = " << w/2.0/Pi << " Hz"<< endl;
    }
    if(D2 > 0)
    {
        double lam = sqrt(D2)/2.0;
            cout << "Aexp(" << lam << "*t)"  << " + Bexp(" << -1*lam << "*t)" << endl;
         eq_fout << "Aexp(" << lam << "*t)"  << " + Bexp(" << -1*lam << "*t)" << endl;
    }
    eq_fout.close();
}

double p1(double x, double y, double z, double pb, double pc, double pd)
{///dy/dx = pb*z + pc*y + pd*x
    return pb*z + pc*y + pd*x;
}

double g1(double x, double y, double z, double gb, double gc, double gd)
{///dz/dx = gb*z + gc*y + gd*x
    return gb*z + gc*y + gd*x;
}

matrix find_incline_mz_points(string Data_Way, string base_file_name)
{
    matrix mz;
        mz.Load_matrix(Data_Way + base_file_name + ".txt");
            cout << mz << endl;

    matrix mz_a = mz;
        mz_a.erase_row(mz.get_height());
            //mz_a.erase_row(mz.get_height());
    double arr_mz_a[mz_a.get_height()-1]; ///на одну точку меньше
    matrix m22;

    for(int i = 0; i < mz_a.get_height(); i++)
    {cout << i << endl;
        m22 = mz.get_row(i);
            m22 = m22.merge_height(mz.get_row(i+1), m22.get_height());
                //m22 = m22.merge_height(mz.get_row(i+2), m22.get_height());
            cout << m22 << endl;
        Find_Incline_Coefficient(m22, arr_mz_a, m22.get_height());
        mz_a.set_data(i,1, *arr_mz_a);
    }
    return mz_a;
}



/*matrix equation_func1(string Data_Way, string file_name, double h, int N, double gb, double gc, double pb, double pc, double x0, double y0, double z0, double shift_angle)///система линейных уравнений
{cout << "equation_func1()" << endl;
//    base_equation(Data_Way, (file_name + base_file_name + "_equation"), gb, gc);
    //matrix mz_a; mz_a.Load_matrix("C:/Users/NiMi/Desktop/Course/7th_stage/Bac/388-sphP/Рассчет/2.268_dmz_a.txt");
    ///
    double k11, k12, k21, k22, k31, k32, k41, k42;
        double x[N]; double y[N]; double z[N];
    ///гр усл /// углы в радианах
    x[0] = x0; y[0] = y0; z[0] = z0;
    for(int i = 0; i < N; i++)
    {        //cout << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;//cout << x[i] << endl;
        ///обновление соответствующим gb, gc от угла
        ///gc = mz_a


        k11 = h*p(x[i], y[i], z[i], pb, pc);                                   k12 = h*g(x[i], y[i], z[i], gb, gc);
        k21 = h*p(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc);       k22 = h*g(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc);
        k31 = h*p(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc);       k32 = h*g(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc);
        k41 = h*p(x[i] + h, y[i] + k31, z[i] + k32, pb, pc);                   k42 = h*g(x[i] + h, y[i] + k31, z[i] + k32, gb, gc);
            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
        x[i+1] = x[i] + h;
    }
        //cout << "eterations finished total iterations = " << N << endl;
    //вообще можно просто проманипулировать указателями
    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N); my = my*180;    ///переход в градусы
            matrix mz; mz.array_to_matrix_column(z, N); mz = mz*180;    ///переход в градусы

    my.column_shift(0, shift_angle);///подвинули полученное решение по углу (балансировочный угол)
        matrix A = mx;
            A = A.merge_width(my, A.get_width());
                A = A.merge_width(mz, A.get_width());
    //A.info();
    //A.Save_matrix(Data_Way + file_name + "num_result.txt"); /// x, z, y <-> x, z, dz
    return A;
};*/

/*matrix equation_func2(string Data_Way, string file_name, double h, int N, double gb0, double gb2, double gc,  double pb, double pc, double x0, double y0, double z0, double shift_angle)///система линейных уравнений
{cout << "equation_func2()" << endl;
//    base_equation(Data_Way, (file_name + base_file_name + "_equation"), gb, gc);
    //matrix mz_a; mz_a.Load_matrix("C:/Users/NiMi/Desktop/Course/7th_stage/Bac/388-sphP/Рассчет/2.268_dmz_a.txt");
    ///
    double gb = gb0;
   // double fr = -0.1;
    double k11, k12, k21, k22, k31, k32, k41, k42;
        double x[N]; double y[N]; double z[N];
    ///гр усл /// углы в радианах
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
    for(int i = 0; i < N; i++)
    {        //cout << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;//cout << x[i] << endl;

        gb = gb0 + gb2*y[i]*y[i];///
            k11 = h*p(x[i], y[i], z[i], pb, pc);                                   k12 = h*g(x[i], y[i], z[i], gb, gc);

        gb = gb0 + gb2*(y[i] + k11/2.0)*(y[i] + k11/2.0);
            k21 = h*p(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc);       k22 = h*g(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc);

        gb = gb0 + gb2*(y[i] + k21/2.0)*(y[i] + k21/2.0);
            k31 = h*p(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc);       k32 = h*g(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc);

        gb = gb0 + gb2*(y[i] + k31/2.0)*(y[i] + k31/2.0);
            k41 = h*p(x[i] + h, y[i] + k31, z[i] + k32, pb, pc);                   k42 = h*g(x[i] + h, y[i] + k31, z[i] + k32, gb, gc);
       // cout << y[i+1]*57.32 << endl;

            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
        x[i+1] = x[i] + h;

     if(fabs(y[i+1])*57.32 >= 90.0)
        {cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }
        //cout << "eterations finished total iterations = " << N << endl;
    //вообще можно просто проманипулировать указателями
    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N); my = my*180/PI;    ///переход в градусы
            matrix mz; mz.array_to_matrix_column(z, N); mz = mz*180/PI;    ///переход в градусы

    my.column_shift(0, shift_angle);///подвинули полученное решение по углу (балансировочный угол)
        matrix A = mx;
            A = A.merge_width(my, A.get_width());
                A = A.merge_width(mz, A.get_width());
    //A.info();
    //A.Save_matrix(Data_Way + file_name + "num_result.txt"); /// x, z, y <-> x, z, dz
    return A;
};*/

/*matrix equation_func3(string Data_Way, string file_name, double h, int N, int fix_max, double gb0, double gb1, double gb2, double gc,  double pb, double pc, double x0, double y0, double z0, double shift_angle)///система линейных уравнений
{cout << "equation_func3()" << endl;
//    base_equation(Data_Way, (file_name + base_file_name + "_equation"), gb, gc);
    //matrix mz_a; mz_a.Load_matrix("C:/Users/NiMi/Desktop/Course/7th_stage/Bac/388-sphP/Рассчет/2.268_dmz_a.txt");
    ///
    double gb = gb0;
   // double fr = -0.1;
    double k11, k12, k21, k22, k31, k32, k41, k42;
        double x[N]; double y[N]; double z[N];
    ///гр усл /// углы в радианах
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
    for(int i = 0; i < N; i++)
    {        //cout << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;//cout << x[i] << endl;

        gb = gb0 + gb1*y[i] + gb2*y[i]*y[i];///
            k11 = h*p(x[i], y[i], z[i], pb, pc);                                   k12 = h*g(x[i], y[i], z[i], gb, gc);

        gb = gb0 + gb1*(y[i] + k11/2.0) + gb2*(y[i] + k11/2.0)*(y[i] + k11/2.0);
            k21 = h*p(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc);       k22 = h*g(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc);

        gb = gb0 + gb1*(y[i] + k21/2.0) + gb2*(y[i] + k21/2.0)*(y[i] + k21/2.0);
            k31 = h*p(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc);       k32 = h*g(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc);

        gb = gb0 + gb1*(y[i] + k31/2.0) + gb2*(y[i] + k31/2.0)*(y[i] + k31/2.0);
            k41 = h*p(x[i] + h, y[i] + k31, z[i] + k32, pb, pc);                   k42 = h*g(x[i] + h, y[i] + k31, z[i] + k32, gb, gc);
       // cout << y[i+1]*57.32 << endl;

            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
        x[i+1] = x[i] + h;

     if(fix_max==1)
        if(fabs(y[i+1])*57.32 >= 90.0)
        {cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }
        //cout << "eterations finished total iterations = " << N << endl;
    //вообще можно просто проманипулировать указателями
    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N); my = my*180/PI;    ///переход в градусы
            matrix mz; mz.array_to_matrix_column(z, N); mz = mz*180/PI;    ///переход в градусы

    my.column_shift(0, shift_angle);///подвинули полученное решение по углу (балансировочный угол)
        matrix A = mx;
            A = A.merge_width(my, A.get_width());
                A = A.merge_width(mz, A.get_width());
    //A.info();
    //A.Save_matrix(Data_Way + file_name + "num_result.txt"); /// x, z, y <-> x, z, dz
    return A;
}*/


/*matrix equation_func4(string Data_Way, string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4, double gc,  double pb, double pc, double x0, double y0, double z0, double shift_angle)///система линейных уравнений
{cout << "equation_func4()" << endl;
//    base_equation(Data_Way, (file_name + base_file_name + "_equation"), gb, gc);
    //matrix mz_a; mz_a.Load_matrix("C:/Users/NiMi/Desktop/Course/7th_stage/Bac/388-sphP/Рассчет/2.268_dmz_a.txt");
    ///
    double gb = gb0;
   // double fr = -0.1;
    double k11, k12, k21, k22, k31, k32, k41, k42;
        double x[N]; double y[N]; double z[N];
    ///гр усл /// углы в радианах
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
    for(int i = 0; i < N; i++)
    {        //cout << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;//cout << x[i] << endl;

        gb = gb0 + gb2*y[i]*y[i] + gb4*pow(y[i],4);///
            k11 = h*p(x[i], y[i], z[i], pb, pc);                                   k12 = h*g(x[i], y[i], z[i], gb, gc);

        gb = gb0 + gb2*(y[i] + k11/2.0)*(y[i] + k11/2.0) + gb4*pow(y[i] + k11/2.0, 4);
            k21 = h*p(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc);       k22 = h*g(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc);

        gb = gb0 + gb2*(y[i] + k21/2.0)*(y[i] + k21/2.0) + gb4*pow(y[i] + k21/2.0, 4);
            k31 = h*p(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc);       k32 = h*g(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc);

        gb = gb0 + gb2*(y[i] + k31/2.0)*(y[i] + k31/2.0) + gb4*pow(y[i] + k31/2.0, 4);
            k41 = h*p(x[i] + h, y[i] + k31, z[i] + k32, pb, pc);                   k42 = h*g(x[i] + h, y[i] + k31, z[i] + k32, gb, gc);
       // cout << y[i+1]*57.32 << endl;

            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
        x[i+1] = x[i] + h;

     if(fix_max==1)
        if(fabs(y[i+1])*57.32 >= 90.0)
        {cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }
        //cout << "eterations finished total iterations = " << N << endl;
    //вообще можно просто проманипулировать указателями
    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N); my = my*180/PI;    ///переход в градусы
            matrix mz; mz.array_to_matrix_column(z, N); mz = mz*180/PI;    ///переход в градусы

    my.column_shift(0, shift_angle);///подвинули полученное решение по углу (балансировочный угол)
        matrix A = mx;
            A = A.merge_width(my, A.get_width());
                A = A.merge_width(mz, A.get_width());
    //A.info();
    //A.Save_matrix(Data_Way + file_name + "num_result.txt"); /// x, z, y <-> x, z, dz
    return A;
}*/

matrix equation_func_VdP_4(string Data_Way, string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4, double gc, double pb, double pc, double x0, double y0, double z0)
{cout << "equation_func_VdP()" << endl;
    double gb = gb0, pd = 0.0, gd = 0.0;
    double k11, k12, k21, k22, k31, k32, k41, k42;
        double x[N]; double y[N]; double z[N];
    ///гр усл /// углы в радианах
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
    cout << gb0 << "\t" << gb2 << "\t" << gb4 << endl;
    for(int i = 0; i < N; i++)
    {        //cout << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;//cout << x[i] << endl;

        gb = gb0 + gb2*y[i]*y[i] + gb4*pow(y[i],4);///
            k11 = h*p1(x[i], y[i], z[i], pb, pc, pd);                                   k12 = h*g1(x[i], y[i], z[i], gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k11/2.0)*(y[i] + k11/2.0) + gb4*pow(y[i] + k11/2.0, 4);
            k21 = h*p1(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc, pd);       k22 = h*g1(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k21/2.0)*(y[i] + k21/2.0) + gb4*pow(y[i] + k21/2.0, 4);
            k31 = h*p1(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc, pd);       k32 = h*g1(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k31/2.0)*(y[i] + k31/2.0) + gb4*pow(y[i] + k31/2.0, 4);
            k41 = h*p1(x[i] + h, y[i] + k31, z[i] + k32, pb, pc, pd);                   k42 = h*g1(x[i] + h, y[i] + k31, z[i] + k32, gb, gc, gd);
       // cout << y[i+1]*57.32 << endl;

            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
        x[i+1] = x[i] + h;

     //if(fix_max==1)
       // if(fabs(y[i+1])*57.32 >= 90.0)
        //{cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }

    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N);
            matrix mz; mz.array_to_matrix_column(z, N);
    matrix A = mx;
        A = A.merge_width(my, A.get_width());
            A = A.merge_width(mz, A.get_width());
    return A;
}

/*matrix equation_func_VdP_4_fr(string Data_Way, string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4, double gc, double pb, double pc, double fr, double x0, double y0, double z0)
{
cout << "equation_func_VdP_fr()" << endl;
    double gb = gb0;
    double k11, k12, k21, k22, k31, k32, k41, k42;
        double x[N]; double y[N]; double z[N];
    ///гр усл /// углы в радианах
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
    cout << gb0 << "\t" << gb2 << "\t" << gb4 << endl;
    for(int i = 0; i < N; i++)
    {        //cout << i << "\t" << x[i] << "\t" << y[i] << "\t" << z[i] << endl;//cout << x[i] << endl;

        gb = gb0 + gb2*y[i]*y[i] + gb4*pow(y[i],4);///
            k11 = h*p(x[i], y[i], z[i], pb, pc);                                   k12 = h*g(x[i], y[i], z[i], gb, gc);

        gb = gb0 + gb2*(y[i] + k11/2.0)*(y[i] + k11/2.0) + gb4*pow(y[i] + k11/2.0, 4);
            k21 = h*p(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc);       k22 = h*g(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc);

        gb = gb0 + gb2*(y[i] + k21/2.0)*(y[i] + k21/2.0) + gb4*pow(y[i] + k21/2.0, 4);
            k31 = h*p(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc);       k32 = h*g(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc);

        gb = gb0 + gb2*(y[i] + k31/2.0)*(y[i] + k31/2.0) + gb4*pow(y[i] + k31/2.0, 4);
            k41 = h*p(x[i] + h, y[i] + k31, z[i] + k32, pb, pc);                   k42 = h*g(x[i] + h, y[i] + k31, z[i] + k32, gb, gc);
       // cout << y[i+1]*57.32 << endl;

            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0 + fr*sign(y[i]);
        x[i+1] = x[i] + h;

     //if(fix_max==1)
       // if(fabs(y[i+1])*57.32 >= 90.0)
        //{cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }

    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N);
            matrix mz; mz.array_to_matrix_column(z, N);
    matrix A = mx;
        A = A.merge_width(my, A.get_width());
            A = A.merge_width(mz, A.get_width());
    return A;
}*/

matrix equation_func_VdP_4_mdyn(string Data_Way, string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4,
                                double gc, double gd, double pb, double pc, double pd, double mdyn, double x0, double y0, double z0)
{
cout << "equation_func_VdP_mdyn()" << endl;
    ///d2y/dx2-mdyn(gb0+gb2*y2+gb4*y4)dy/dx-x
        ///dz/dx = gb*z+gc*y + gd*x
        ///dy/dx = pb*z+pc*y + pd*x /// p(x,y,z, pb, pc)
    ///d2y/dx2-mdyn(1+e2*y2+e4*y4)dy/dx+x
    double gb = gb0;
    double k11, k12, k21, k22, k31, k32, k41, k42;
        double *x = new double [N+1];
        double *y = new double [N+1];
        double *z = new double [N+1];

        double f = 0.0, g = 0.0;
        double *lam1Re = new double [N+1];
        double *lam1Im = new double [N+1];
        double *lam2Re = new double [N+1];
        double *lam2Im = new double [N+1];
    ///гр усл /// углы в радианах
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
            cout << "dz/dx = "<< mdyn <<  "*(" << gb4*mdyn << "*y^4 + " << gb2*mdyn << "*y^2 + "  << gb0*mdyn << ")*z + " << gc << "*y = 0" << endl;
            cout << "dy/dx = " << pb << "*z + " << pc << "*y = 0" << endl;
    //cout << gb0 << "\t" << gb2 << "\t" << gb4 << endl;
    for(int i = 0; i < N; i++)
    {//cout << i<< endl;
        gb = gb0 + gb2*y[i]*y[i] + gb4*pow(y[i],4); gb *= mdyn;///
            k11 = h*p1(x[i], y[i], z[i], pb, pc, pd);                                   k12 = h*g1(x[i], y[i], z[i], gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k11/2.0)*(y[i] + k11/2.0) + gb4*pow(y[i] + k11/2.0, 4); gb *= mdyn;
            k21 = h*p1(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc, pd);       k22 = h*g1(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k21/2.0)*(y[i] + k21/2.0) + gb4*pow(y[i] + k21/2.0, 4); gb *= mdyn;
            k31 = h*p1(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc, pd);       k32 = h*g1(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k31/2.0)*(y[i] + k31/2.0) + gb4*pow(y[i] + k31/2.0, 4); gb *= mdyn;
            k41 = h*p1(x[i] + h, y[i] + k31, z[i] + k32, pb, pc, pd);                   k42 = h*g1(x[i] + h, y[i] + k31, z[i] + k32, gb, gc, gd);

            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
        x[i+1] = x[i] + h;

        ///анализ собственных чисел матрицы якоби lam1(Re, Im) и lam2(Re, Im)
        f = gb0 + gb2*y[i]*y[i] + gb4*pow(y[i],4);
        g = gc;
        equation_2order(1, -f, -g, lam1Re[i], lam1Im[i], lam2Re[i], lam2Im[i]);

        if(fabs(y[i+1])>= 3.14*10)
        {cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }
    cout << "MST = " << gc << "\tf0 = " << sqrt(-gc)/PI/2.0 << "\tw0 = " << sqrt(-gc) << endl;

    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N);
            matrix mz; mz.array_to_matrix_column(z, N);

    matrix l1Re; l1Re.array_to_matrix_column(lam1Re, N);
        matrix l1Im; l1Im.array_to_matrix_column(lam1Im, N);
            matrix l2Re; l2Re.array_to_matrix_column(lam2Re, N);
                matrix l2Im; l2Im.array_to_matrix_column(lam2Im, N);

    matrix A = mx;
        A = A.merge_width(my);
            A = A.merge_width(mz);
    A = A.merge_width(l1Re);
        A = A.merge_width(l1Im);
            A = A.merge_width(l2Re);
                A = A.merge_width(l2Im);
    delete [] x;    delete [] y;    delete [] z;

    delete [] lam1Re;    delete [] lam1Im;    delete [] lam2Re;    delete [] lam2Im;
    return A;
}

matrix equation_func_VdP_4_mdyn_poly_mst(string Data_Way, string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4,
                                          double gc0, double gc2, double gc4, double gd, double pb, double pc, double pd,
                                            double mdyn, double x0, double y0, double z0)
{
cout << "equation_func_VdP_mdyn_poly_mst()" << endl;
    ///d2y/dx2-mdyn(gb0+gb2*y2+gb4*y4)dy/dx-x
        ///dz/dx = gb*z+gc*y + gd*x
        ///dy/dx = pb*z+pc*y + pd*x /// p(x,y,z, pb, pc)
    ///d2y/dx2-mdyn(1+e2*y2+e4*y4)dy/dx+x
    double gb = gb0, gc = gc0;
    double k11, k12, k21, k22, k31, k32, k41, k42;

    ///решение
        double* x = new double[N]; double* y = new double [N]; double* z = new double[N];
    double f = 0.0, g = 0.0;
        double *lam1Re = new double [N];double *lam1Im = new double [N];
        double *lam2Re = new double [N];double *lam2Im = new double [N];

    ///гр усл
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
            cout << "dz/dx = "<< mdyn <<  "*(" << gb4*mdyn << "*y^4 + " << gb2*mdyn << "*y^2 + "  << gb0*mdyn << ")*z + "
                << gc0 << "*y + " << gc2 << "*y^3 + " << gc4 << "*y^5 = 0" << endl;
            cout << "dy/dx = " << pb << "*z + " << pc << "*y = 0" << endl;
    //cout << gb0 << "\t" << gb2 << "\t" << gb4 << endl;
    for(int i = 0; i < N; i++)
    {
        gb = gb0 + gb2*y[i]*y[i] + gb4*pow(y[i],4); gb *= mdyn;///
        gc = gc0 + gc2*y[i]*y[i] + gc4*pow(y[i],4);///исправлено на gc4 вместо gb4
            k11 = h*p1(x[i], y[i], z[i], pb, pc, pd);                                   k12 = h*g1(x[i], y[i], z[i], gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k11/2.0)*(y[i] + k11/2.0) + gb4*pow(y[i] + k11/2.0, 4); gb *= mdyn;
        gc = gc0 + gc2*(y[i] + k11/2.0)*(y[i] + k11/2.0) + gc4*pow(y[i] + k11/2.0, 4);
            k21 = h*p1(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc, pd);       k22 = h*g1(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k21/2.0)*(y[i] + k21/2.0) + gb4*pow(y[i] + k21/2.0, 4); gb *= mdyn;
        gc = gc0 + gc2*(y[i] + k21/2.0)*(y[i] + k21/2.0) + gc4*pow(y[i] + k21/2.0, 4);
            k31 = h*p1(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc, pd);       k32 = h*g1(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc, gd);

        gb = gb0 + gb2*(y[i] + k31/2.0)*(y[i] + k31/2.0) + gb4*pow(y[i] + k31/2.0, 4); gb *= mdyn;
        gc = gc0 + gc2*(y[i] + k31/2.0)*(y[i] + k31/2.0) + gc4*pow(y[i] + k31/2.0, 4);
            k41 = h*p1(x[i] + h, y[i] + k31, z[i] + k32, pb, pc, pd);                   k42 = h*g1(x[i] + h, y[i] + k31, z[i] + k32, gb, gc, gd);

        y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;           z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
        x[i+1] = x[i] + h;

        ///анализ собственных чисел матрицы якоби lam1(Re, Im) и lam2(Re, Im)
        f = gb0 + gb2*y[i]*y[i] + gb4*pow(y[i],4);
        g = gc0 + gc2*y[i]*y[i] + gc4*pow(y[i],4);
            equation_2order(1.0, -f, -g, lam1Re[i], lam1Im[i], lam2Re[i], lam2Im[i]);
        /*///+D
        if(f*f+4.0*g >= 0.0)
        {
            lam1Re[i] = (f+sqrt(f*f+4.0*g))/2.0;
            lam1Im[i] = 0.0;
        }else{
            lam1Re[i] = f/2.0;
            lam1Im[i] = sqrt(-1.0*(f*f+4.0*g))/2.0;;
        }
        ///-D
        if(f*f+4.0*g >= 0.0)
        {
            lam2Re[i] = (f-sqrt(f*f+4.0*g))/2.0;
            lam2Im[i] = 0.0;
        }else{
            lam2Re[i] = f/2.0;
            lam2Im[i] = -1.0*sqrt(-1.0*(f*f+4.0*g))/2.0;;
        }*/

     //if(fix_max==1)
        if(fabs(y[i+1])>= 3.14*10)///сомнительно
        {cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }

    matrix mx; mx.array_to_matrix_column(x, N);
        matrix my; my.array_to_matrix_column(y, N);
            matrix mz; mz.array_to_matrix_column(z, N);

    matrix l1Re; l1Re.array_to_matrix_column(lam1Re, N);
        matrix l1Im; l1Im.array_to_matrix_column(lam1Im, N);
            matrix l2Re; l2Re.array_to_matrix_column(lam2Re, N);
                matrix l2Im; l2Im.array_to_matrix_column(lam2Im, N);

    matrix A = mx;
        A = A.merge_width(my, A.get_width());
            A = A.merge_width(mz, A.get_width());

    A = A.merge_width(l1Re);
        A = A.merge_width(l1Im);
            A = A.merge_width(l2Re);
                A = A.merge_width(l2Im);

    delete [] x; delete [] y; delete [] z;

    delete [] lam1Re;    delete [] lam1Im;    delete [] lam2Re;    delete [] lam2Im;
    return A;
}

matrix equation_func_polynomial4order_f_g(string Data_Way, string file_name, double h, int N, int fix_max,
                                          function<double(double fa0, double fa2, double fa4, double f_arg)> f_poly,
                                          function<double(double ga0, double ga2, double ga4, double g_arg)> g_poly,
                                          double gb0, double gb2, double gb4,double gc0, double gc2, double gc4,
                                          double gd, double pb, double pc, double pd,double mdyn, double x0, double y0, double z0)
{
cout << "equation_func_VdP_mdyn_poly_mst()" << endl;
    ///d2y/dx2-mdyn(gb0+gb2*y2+gb4*y4)dy/dx-x ///dz/dx = gb*z+gc*y + gd*x///dy/dx = pb*z+pc*y + pd*x /// p(x,y,z, pb, pc)
    double gb = gb0, gc = gc0;
    double k11, k12, k21, k22, k31, k32, k41, k42;
    ///решение
        double* x = new double[N]; double* y = new double [N]; double* z = new double[N];
    double f = 0.0, g = 0.0;
        double *lam1Re = new double [N];double *lam1Im = new double [N];
        double *lam2Re = new double [N];double *lam2Im = new double [N];
    ///гр усл
    x[0] = x0; y[0] = y0; z[0] = z0;
    int out_index = 0;
        cout << "equation:\n\tdz/dx = "<< mdyn <<  "*(" << gb4*mdyn << "*y^4 + " << gb2*mdyn << "*y^2 + "  << gb0*mdyn << ")*z + "
            << gc0 << "*y + " << gc2 << "*y^3 + " << gc4 << "*y^5 = 0" << endl;
        cout << "\tdy/dx = " << pb << "*z + " << pc << "*y = 0" << endl;
    ///dy/dx = pb*z + pc*y + pd*x///dz/dx = gb*z + gc*y + gd*x
    function<double(double,double,double, double,double,double)> Eq_func = [](double argX, double argY, double argZ, double coefZ, double coefY, double coefX)
        {return coefZ*argZ + coefY*argY + coefX*argX;};

    for(int i = 0; i < N; i++)
    {
        k11 = h*Eq_func(x[i], y[i], z[i], pb, pc, pd);
            gb = f_poly(gb0, gb2, gb4, y[i])*mdyn;
            gc = g_poly(gc0, gc2, gc4, y[i]);
        k12 = h*Eq_func(x[i], y[i], z[i], gb, gc, gd);

        k21 = h*Eq_func(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, pb, pc, pd);
            gb = f_poly(gb0, gb2, gb4, y[i] + k11/2.0)*mdyn;
            gc = g_poly(gc0, gc2, gc4, y[i] + k11/2.0);
        k22 = h*Eq_func(x[i] + h/2.0, y[i] + k11/2.0, z[i] + k12/2.0, gb, gc, gd);

        k31 = h*Eq_func(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, pb, pc, pd);
            gb = f_poly(gb0, gb2, gb4, y[i] + k21/2.0)*mdyn;
            gc = g_poly(gc0, gc2, gc4, y[i] + k21/2.0);
        k32 = h*Eq_func(x[i] + h/2.0, y[i] + k21/2.0, z[i] + k22/2.0, gb, gc, gd);

        k41 = h*Eq_func(x[i] + h, y[i] + k31, z[i] + k32, pb, pc, pd);
            gb = f_poly(gb0, gb2, gb4, y[i] + k31/2.0)*mdyn;
            gc = g_poly(gc0, gc2, gc4, y[i] + k31/2.0);
        k42 = h*Eq_func(x[i] + h, y[i] + k31, z[i] + k32, gb, gc, gd);

        z[i+1] = z[i] + (k12 + 2.0*k22 + 2.0*k32 + k42)/6.0;
            y[i+1] = y[i] + (k11 + 2.0*k21 + 2.0*k31 + k41)/6.0;
                x[i+1] = x[i] + h;

        ///анализ собственных чисел матрицы якоби lam1(Re, Im) и lam2(Re, Im)
        f = f_poly(gb0, gb2, gb4, y[i]);///тут можно просто было взять gb и gc
        g = g_poly(gc0, gc2, gc4, y[i]);
            equation_2order(1.0, -f, -g, lam1Re[i], lam1Im[i], lam2Re[i], lam2Im[i]);

        if(fabs(y[i+1])>= 3.14*10)///результат уходит за реальные значения - остановка счета
            {cout << y[i+1]*57.32 << endl; y[i+1] = 0.0; out_index = N; N=i+1; i = out_index; break;}
    }

    matrix result, buff_column;

    buff_column = buff_column.matrix_column_from_array(x, N);

    result = buff_column;
        result = result.merge_width(buff_column.matrix_column_from_array(y,N));
            result = result.merge_width(buff_column.matrix_column_from_array(z,N));
    result = result.merge_width(buff_column.matrix_column_from_array(lam1Re, N));
        result = result.merge_width(buff_column.matrix_column_from_array(lam1Im, N));
            result = result.merge_width(buff_column.matrix_column_from_array(lam2Re, N));
                result = result.merge_width(buff_column.matrix_column_from_array(lam2Im, N));

    delete [] x; delete [] y; delete [] z;
        delete [] lam1Re;    delete [] lam1Im;    delete [] lam2Re;    delete [] lam2Im;
    return result;
}

matrix equation_func_implicit_Milan(string Data_Way, string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4,
                                          double gc0, double gc2, double gc4, double gd, double pb, double pc, double pd,
                                            double mdyn, double x0, double y0, double z0)
{
    ///implicit milan schime
    matrix m;
    return m;
}

