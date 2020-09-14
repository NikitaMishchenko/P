#include <iostream>
#include <cmath>
#include <string>

#include "manager.h"
#include "equation_num.h"

using namespace std;
/// DIFFERENCIAL EQUATIONS

///IMPLICIT
void Calc_VdP_equation_4(string Data_Way, string base_file_name, int fix_max)
{
    int N_steps;
        ///¬»ƒ ”–ј¬Ќ≈Ќ»я
        double T, h, lam, mu;
        double gb, gc, pb, pc, x0, y0, z0, shift_angle = 0, shift_time = 0;
        matrix C_RES, RES, d_RES, Env_RES, Ampl_RES(1,4), Empt;///Ampl_RES //lam, mu, rho-, rho+

        string s_buff = "", new_name;
        ifstream fin(Data_Way + "VdP_4_param.txt");

    while(new_name != "START")
        fin >> new_name;///buff
    fin >> new_name; ///data
    cout << "new_name = " << new_name << endl;
    while(new_name != "END")
    {
            fin >>
                T >> h >>
                    x0 >> y0 >> z0 >> ///начальные значени€ в градусах!
                            lam >> mu;/// коэффициенты уравнени€ при производной пор€дка 1, 1 и (x*x), 0 (при 2 = 1)
            N_steps = T/h;///количество шагов по времени и указанному шагу

        ///d2y/dx2-(lam+mu*y2-y4)dy/dx+x
        ///d2y/dx2+(-lam-mu*y2+y4)dy/dx+x
        pb = 1.0; pc = 0.0;
        C_RES = equation_func_VdP_4(Data_Way, base_file_name, h, N_steps, fix_max, -lam, -mu, 1.0, 1, pb, pc, x0, y0, z0);;
            C_RES.Save_matrix(Data_Way + new_name + "_VdP_eq.txt");

        RES = C_RES.get_column(0);
            RES = RES.merge_width(C_RES.get_column(1), RES.get_width());
        d_RES = C_RES.get_column(0);
            d_RES = d_RES.merge_width(C_RES.get_column(2), RES.get_width());

        Env_RES = Get_Envelope_top(RES, d_RES);
            Env_RES.Save_matrix(Data_Way + new_name + "_envelop_eq.txt");
        ///амплитуда минимальна€ и максимальна€
        Ampl_RES.set_data(0,0, lam);
        Ampl_RES.set_data(0,1, mu);
        Ampl_RES.set_data(0,2, Env_RES.get_min_el_column(1));
        Ampl_RES.set_data(0,3, Env_RES.get_max_el_column(1));
            Ampl_RES.Save_matrix_end(Data_Way + new_name + "_ampl_eq.txt");

        fin >> new_name ;/// следующий или конец?
        C_RES = Empt;
    }
    fin.close();
}

void Calc_VdP_equation_4_mdyn(string Data_Way, string base_file_name, string mode, int fix_max)
{cout << "Calc_VdP_equation_4_mdyn(string Data_Way, string base_file_name, string mode, int fix_max)\n";
    int N_steps;
        ///¬»ƒ ”–ј¬Ќ≈Ќ»я
        double T, h, e0, e2, e4, mdyn, gc0, gc2, gc4; ///врем€ счета, шаг счета, коэффициенты.
        double gb, gc, gd, pb, pc, pd, x0, y0, z0, shift_angle = 0, shift_time = 0;
        matrix C_RES, RES, d_RES, Env_RES, Ampl_RES(1,4), const_column;///Ampl_RES //lam, mu, rho-, rho+

        string s_buff = "", new_name;
        ifstream fin(Data_Way + "VdP_4_mdyn_param.txt");

    while(new_name != "START")
        fin >> new_name;///buff
    fin >> new_name; ///data
        cout << "DATA = " << new_name << endl;
    while(new_name != "END")
    {
            fin >>  T >> h >>
                    x0 >> y0 >> z0 >> ///начальные значени€ в градусах!
                            e0 >> e2 >> e4 >> mdyn >>
                                gc0 >> gc2 >> gc4;/// коэффициенты уравнени€ при производной пор€дка 1, 1 и (x*x), 0 (при 2 = 1)
            N_steps = T/h;///количество шагов по времени и указанному шагу

        gb = 0.0; gc = -1.0; gd = 0.0; ///dz/dx = gb()*z - x
        pb = 1.0; pc = 0.0; pd = 0.0; ///dy/dx = z

        if (mode == "0")
            C_RES = equation_func_VdP_4_mdyn(Data_Way, base_file_name, h, N_steps, fix_max, e0, 1.0*e2, 1.0*e4, gc, gd, pb, pc, pd, mdyn, x0, y0, z0);
        if (mode == "exp_eq")
            C_RES = equation_func_VdP_4_mdyn(Data_Way, base_file_name, h, N_steps, fix_max, e0, 1.0*e2, 1.0*e4, gc0, gd, pb, pc, pd, mdyn, x0, y0, z0);
        if (mode == "VdP4_mdyn_poly_mst")
            C_RES = equation_func_VdP_4_mdyn_poly_mst(Data_Way, base_file_name, h, N_steps, fix_max,e0, 3.0*e2, 5.0*e4,gc0, gc2, gc4, gd, pb, pc, pd, mdyn,x0, y0, z0);
        if (mode == "advanced_RK4")
        {
            cout << "advanced_RK4\n\tf = Poly[4], g = Poly[1]\n";
                C_RES = equation_func_polynomial4order_f_g(Data_Way, base_file_name, h, N_steps, fix_max,
                                                   [](double fa0, double fa2, double fa4, double f_arg){return fa0 + fa2*pow(f_arg,2) + fa4*pow(f_arg,4);},
                                                   [](double ga0, double ga2, double ga4, double g_arg){return ga0 + ga2*pow(g_arg,2) + ga4*pow(g_arg,4);},
                                                   e0, 3.0*e2, 5.0*e4,gc0, gc2, gc4,
                                                   gd, pb, pc, pd, mdyn, x0, y0, z0);
        }

        ///сохранение результатов численного решени€ уравнени€ и сч матрицы якоби: x,y,z,lam1Re,lam1Im,lam2Re,lam2Im
        C_RES.Save_matrix(Data_Way + new_name + "_VdP_eq" + ".txt");

        ///сохранение в форме oscillation данных дл€ обработки: count1, flow3, angle1, dangle1, freq1, time1 = 8
        const_column.make_matrix_data(C_RES.get_height(),1,1.0);/// имеетс€ врем€, угол, скорость
        RES = RES.merge_width(const_column);///count
        RES = RES.merge_width(const_column);///flow1
        RES = RES.merge_width(const_column);///flow2
        RES = RES.merge_width(const_column);///flow3
        RES = RES.merge_width(C_RES.get_column(1));///angle
        RES = RES.merge_width(C_RES.get_column(2));///dangle
        RES = RES.merge_width(const_column);///freq
        RES = RES.merge_width(C_RES.get_column(0));///time
        RES.Save_matrix(Data_Way + new_name + "_r" + ".txt");

        C_RES.delete_data();///освобождение дл€ следующей задачи
        fin >> new_name ;/// следующий или конец?
    }
    fin.close();
}
