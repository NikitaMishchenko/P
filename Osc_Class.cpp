#include <iostream>

#include <cmath>

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <string>

#include "matrix.h"
#include "Approximation.h"

#include "Osc_Class.h"

using namespace std;

const double PI = 3.14159265359;


oscillation::oscillation()
{
    //matrix t(0,1);
      //  time = t;
    matrix a(0,1);
        angle = a;
    matrix da(0,1);
        dangle = da;
    matrix fl(0, 3);
        flow = fl;

    matrix fr(0,1);
        freq = fr;

    matrix b_a(0,1);
        balance_angle = b_a;

    calibr_angle = 0;
    data_source = "no_source";
    model_name = "";
    I = 0;    S = 0;    l = 0;
}

oscillation::oscillation(int h)
{
    //matrix t(h, 1);
   //     time = t;
    matrix a(h, 2);
        angle = a;
    matrix da(h, 1);
        dangle = da;
    matrix fl(h, 3);
        flow = fl;
    matrix fr(h, 1);
        freq = fr;
    matrix b_a(h,1);
        balance_angle = b_a;

    calibr_angle = 0;
    data_source = "no_source";
    model_name = "";
    I = 0;    S = 0;    l = 0;
};

oscillation::oscillation(string file_name)
{
    matrix buff;
        buff.Load_matrix(file_name);
    if(buff.get_width() >= 6)
    {
        angle = buff.get_column(0);
        angle = angle.merge_width(buff.get_column(1), angle.get_width());
        dangle = buff.get_column(2);
        flow = buff.get_column(3);
            flow = flow.merge_width( buff.get_column(4), flow.get_width());
                flow = flow.merge_width( buff.get_column(5), flow.get_width());
        calibr_angle = 0;
        data_source = file_name;
    }else{
        oscillation O;
        *this = O;
    }
    model_name = "";
    data_source = "";
    matrix b_a(angle.get_height(),1);
        balance_angle = b_a;
    I = 0;    S = 0;    l = 0;
};

double oscillation::get_init_angle()
    {return init_angle;}
//void oscillation::set_balance_angle(double a)
 //   { balance_angle = a;}
//double oscillation::get_balance_angle()
//    { return balance_angle;}

///main
double oscillation::calc_balance_angle()
{
    double b_angle = 0.0;
    /*angle.info();
    for(int i = 0; i < angle.get_height(); i++)
    {
        b_angle += angle.get_data(i,0);
    }
    angle = angle/angle.get_height();*/
    b_angle = angle.get_data(angle.get_height()-1, 0);
    return b_angle;
}

/*matrix oscillation::count_freq()
{
    matrix W(angle.get_height(), 1);
        double h = time.get_data(1, 0) -  angle.get_data(0, 0);

    W.set_data(0, 0, (-3*angle.get_data(0, 0) + 4*angle.get_data(1, 0) - angle.get_data(2, 0))/2/h );
    for(int i = 1; i < angle.get_height()-1; i++)
        W.set_data(i, 0, (angle.get_data(i+1, 0) - angle.get_data(i-1, 0))/2/h );
    W.set_data(W.get_height()-1, 0, (angle.get_data(A.get_height()-3, 0) - 4*angle.get_data(W.get_height()-2, 0) + 3*angle.get_data(W.get_height()-1, 0))/2/h );

    return W;
}*/

void oscillation::angle_shift_balance(int index_column)
{cout << "oscillation::angle_shift_balance()\n";
    for(int i = 0; i < balance_angle.get_height(); i++)
        {
        //cout << angle.get_data(i,0) << "\t" << this->get_balance_angle_data(i,index_column) << "\t";
            angle.set_data(i,0, angle.get_data(i,0) - this->get_balance_angle_data(i,index_column));
            //cout << angle.get_data(i,0) << endl;
        }
}

void oscillation::angle_shift_balance()
    {this->angle_shift_balance(0);}

matrix oscillation::count_freq(matrix t, matrix A)
{//cout << "oscillation::count_freq(matrix t, matrix A)\n";
    matrix W;//matrix W(A.get_height(), 1);
    //double h = t.get_data(1, 0) -  t.get_data(0, 0);
    W = derevative_3_dot(A, t);
    //cout << "freq";

    /*W.set_data(0, 0, (-3*A.get_data(0, 0) + 4*A.get_data(1,0) - A.get_data(2, 0))/2/h);
    for(int i = 1; i < A.get_height()-1; i++)
        W.set_data(i, 0, (A.get_data(i+1, 0) - A.get_data(i-1, 0))/2/h );
    W.set_data(W.get_height()-1, 0, (A.get_data(A.get_height()-3, 0) - 4*A.get_data(W.get_height()-2, 0) + 3*A.get_data(W.get_height()-1,0))/2/h );*/
    return W.get_column(1);
};


void oscillation::callibrate_to_zero()
{
    if(calibr_angle != 0)
        angle.column_shift(0, calibr_angle);
}

void oscillation::callibr_angle()
{
    double sum;
    for(int i = 0; i < angle.get_height(); i++)
        sum += angle.get_data(i,0);
    sum = sum/angle.get_height();
    for(int i = 0; i < angle.get_height(); i++)
        angle.set_data(i,0, angle.get_data(i,0) - sum);
}

void oscillation::shift_time(double time_shift) ///flow тоже нужно сдвигать
    {angle.column_shift(0, time_shift);};

void oscillation::shift_angle(double angle_shift)
    {angle.column_shift(0, angle_shift);};

void oscillation::callibrate_to_radian()
{
    angle = angle*(PI/180.0);
    calibr_angle *=(PI/180.0);
    dangle = dangle*(PI/180.0);
};

void oscillation::callibrate_to_degree()
{
    angle = angle*(180.0/PI);
    calibr_angle *= (180.0/PI);
    dangle = dangle*(180.0/PI);
};


///primary
void oscillation::set_time(matrix A)
    {time = A;};
matrix oscillation::get_time()
    {return time;};


void oscillation::set_angle(matrix A)
    {angle = A;};
void oscillation::set_angle_data(int i, int j, double d)
    {angle.set_data(i,j,d);};
matrix oscillation::get_angle()
    {return angle;};

double oscillation::get_calibr_angle()
    {return calibr_angle;};
void oscillation::set_calibr_angle(double a)
    {calibr_angle = a;};

void oscillation::set_dangle(matrix A)
    {dangle = A;};
matrix oscillation::get_dangle()
    {return dangle;};

void oscillation::set_freq(matrix f)
    {freq = f;}
matrix oscillation::get_freq()
    {return freq;};

void oscillation::set_flow(matrix f)
    {flow = f;};
matrix  oscillation::get_flow()
    {return flow;};
double oscillation::get_flow_Re(int i)///Reynolds num in 10^6
{
    if(flow.get_width() > 3)
        return flow.get_data(i, 3);

return flow.get_data(i,0);
}
double oscillation::get_flow_M(int i)///Mach num
{
    if(flow.get_width() > 3)
        return flow.get_data(i, 4);
return flow.get_data(i,1);
}
double oscillation::get_flow_v(int i)///veloc
{
    if(flow.get_width() > 3)
        return flow.get_data(i, 4);
return 0;
}
double oscillation::get_flow_q(int i)///head in Pa
{
        if(flow.get_width() >3)
            return flow.get_data(i, 6);
return flow.get_data(i,2);
}


void oscillation::set_flow_data(int i, int j, double data)
    {flow.set_data(i,j,data);};
double oscillation::get_flow_data(int i, int j)
    {return flow.get_data(i,j);};

void oscillation::set_balance_angle(matrix b_a)
    {  balance_angle = b_a; }
void oscillation::set_balance_angle_data(int i, int j, double d)
    {   balance_angle.set_data(i,j,d);};
void oscillation::set_balance_angle_data(int i, double d)
    { balance_angle.set_data(i, 0, d);}
matrix oscillation::get_balance_angle()
    {return balance_angle;}
double oscillation::get_balance_angle_data(int i, int j)
    {return balance_angle.get_data(i,j);}
double oscillation::get_balance_angle_data(int i)
    {return balance_angle.get_data(i,0);}


void oscillation::set_data_source(string data)
    {data_source = data;};
string oscillation::get_data_source()
    {return data_source;}


void oscillation::set_model(string n_mn)
    {model_name = n_mn;};
string oscillation::get_model()
    {return model_name;}
void oscillation::set_model_I(double d)
    {I = d;}
double oscillation::get_model_I()
    {return I;};
void oscillation::set_model_S(double d)
    {S = d;}
double oscillation::get_model_S()
    {return S;};
void oscillation::set_model_l(double d)
    {l = d;}
double oscillation::get_model_l()
    {return l;};

void oscillation::cut_from_to(double time1, double time2)///сдвинул начальную точку влево на 0.01
{cout << "oscillation::cut_from_to(double time1, double time2)\t";
    //cout << time1 <<"\t to " << time2 << endl;
    time1 -=0.01;
    double t = 0.0;
    if(time1 > time2)
        {t = time2; time2 = time1; time1 = t;}

    ///если указали больше или меньше будут даны имеющегося размера
    if(time2 > time.get_data(time.get_height()-1, 0))
        time2 = time.get_data(time.get_height()-1, 0);
    if(time1 < time.get_data(0, 0))
        time1 = time.get_data(0, 0);

    ///находим индексы границ для указанных времен
    int h1 = 0, h2 = time.get_height()-1;// cout << "\t\t" << time1 << "\t" << time2 << "\nFORRRRRRRRRRRRRRRRRRRRRRRRRRRRR" << endl; cin >> h1;
    for(int i = 0; i < time.get_height()-1; i++)
    {//cout << i << "\t" << h1 << "\t" << h2 << "\t" << time.get_data(i,0) << endl;
        if(time.get_data(i,0) <= time1)
            h1 = i;
        if(time.get_data(i,0) <= time2)
        {
            h2 = i;
        }//else{i = time.get_height()-1;}
    }
    //cout << h1 << "\t" << h2 << endl;
    //cin >> h1;
    time = time.get_matrix_part(h1,h2,0,1);//cout << "time\n";
    angle = angle.get_matrix_part(h1,h2,0,1);//cout << "angle\n";
    dangle = dangle.get_matrix_part(h1,h2,0,1);//cout << "dangle\n";
    freq = freq.get_matrix_part(h1,h2,0,1); //cout << "freq\n";
    flow = flow.get_matrix_part(h1,h2,0,3);//cout << "flow\n";
};

oscillation oscillation::oscillation_cut_from_to(double time1, double time2)
{cout << "oscillation::cut_from_to(double time1, double time2)\t";
    cout << time1 << time2 << endl;
    oscillation A;
    double t = 0.0;
    if(time1 > time2)
        {t = time2; time2 = time1; time1 = t;}
    //angle.info();
    //cout << angle.get_data(0, 0) << "\t" << angle.get_data(angle.get_height()-1, 0) << "\t" << angle.get_data(0, 1) << "\t" << angle.get_data(angle.get_height()-1, 1) << endl;
    ///если указали больше или меньше будут даны имеющегося размера
    if(time2 > time.get_data(time.get_height()-1, 0))
        time2 = time.get_data(time.get_height()-1, 0);
    if(time1 < time.get_data(0, 0))
        time1 = time.get_data(0, 0);

    ///находим индексы границ для указанных времен
    int h1 = 0, h2 = 0; //cout << time1 << "\t" << time2 << endl;
    for(int i = 0; i < time.get_height()-1; i++)
    {//cout << i << "\t" << h1 << "\t" << h2 << endl;
        if(time.get_data(i,0) <= time1)
            h1 = i;
        if(time.get_data(i,0) <= time2)
        {
            h2 = i;
        }else{i = time.get_height()-1;}
    }//cout << h1 << h2 << endl;
    A.time = time.get_matrix_part(h1,h2,0,1);//cout << "time\n";
    A.angle = angle.get_matrix_part(h1,h2,0,1);//cout << "angle\n";
    A.dangle = dangle.get_matrix_part(h1,h2,0,1);//cout << "dangle\n";
    A.freq = freq.get_matrix_part(h1,h2,0,1); //cout << "freq\n";
    A.flow = flow.get_matrix_part(h1,h2,0,3);//cout << "flow\n";
    return A;
};


matrix oscillation::w_harm()
{
    matrix dA = derevative_3_dot(this->get_angle(),1,0);
    int arr[dA.get_height()];
        int j = 0;
            for(int i = 1; i < dA.get_height()-1; i++)
                if(dA.get_data(i-1,1) <= 0.0)
                    if(dA.get_data(i,1) >= 0.0)
                    {
                        arr[j] = i;
                        j++;
                    }
    matrix E(j-1, 2);
    for(int i = 0; i < E.get_height(); i++)
    {
        E.set_data(i, 0, (dA.get_data(arr[i+1], 0) + dA.get_data(arr[i], 0))/2.0);
        E.set_data(i, 1, 1.0/(dA.get_data(arr[i+1], 0) - dA.get_data(arr[i], 0)));
    }
    return E;
};

matrix oscillation::w_harm1()
{//cout << "oscillation::w_harm1()" << endl;
    matrix dA = dangle.merge_width(time,0);
        //this->info();
    int* arr = new int[dA.get_height()]; /// массив-карта
        int j = 0;
            for(int i = 1; (j < dA.get_height()-1)  && (i < dA.get_height()-1); i++)///j instead i
            {//cout << i << "\t" << j <<endl;
                if(dA.get_data(i-1,1) <= 0.0)
                    if(dA.get_data(i,1) >= 0.0)
                    {
                        arr[j] = i;
                        j++;
                    }
                if(dA.get_data(i-1,1) >= 0.0)
                    if(dA.get_data(i,1) <= 0.0)
                    {
                        arr[j] = i;
                        j++;
                    }
            }

    matrix E(j-1, 2);
    for(int i = 0; i < E.get_height(); i++)
    {//cout << i << "\t" << j <<endl;

        E.set_data(i, 0, dA.get_data(arr[i], 0)/*(dA.get_data(arr[i+1], 0) + dA.get_data(arr[i], 0))/2.0*/);
        E.set_data(i, 1, 1.0/(dA.get_data(arr[i+1], 0) - dA.get_data(arr[i], 0))/2.0);
    }//E.info();

    matrix E1(this->angle.get_height(),2);

    for(int i = 0, j = 0; j < E.get_height(),i < angle.get_height()-1; )
    {///cout << i << "\t" << angle.get_data(i+1,0) << "\t" << E.get_data(j,0) << endl;
        while(time.get_data(i,0) <= E.get_data(j,0) && i <= angle.get_height()-1)
            {E1.set_data(i,1, E.get_data(j,1));i++;}
               j++;
        if(j == E.get_height()-1)
            while (i<angle.get_height())
            {E1.set_data(i,1, E.get_data(j,1));i++;}
    }
    ///last
    //E1.info();
    delete [] arr;
    return E1;
}
void oscillation::middle_angle_timeavg_halfT(const matrix& env_top,const matrix& env_bot, matrix& middle_angle)
{
    matrix result,res_str(1,2);
    int i = 0, ///счетчик углов
        j1 = 0, ///счетчик огибающих top
        j2 = 0, ///счетчик огибающих bot
        sum_counts = 0;///число отсчетов угла в периоде
    double sum_half = 0.0;

    while( j1 < env_top.get_height()-1 &&  j2 < env_bot.get_height()-1)
    {//cout << j << "\t" << env_bot.get_height()-1 << "\t" << env_top.get_data(j, 0) << "\t" << env_bot.get_data(j, 0);
        if(env_top.get_data(j1, 0) < env_bot.get_data(j2, 0))///верхний раньше нижнего
        {
            while(time.get_data(i,0) < env_top.get_data(j1, 0))
                i++;
            ///суммируем углы между точками огибающей
            while(time.get_data(i,0) < env_bot.get_data(j2, 0))
            {//cout << "\t" << time.get_data(i,0) << "\t" << j << "\t" << i;
                sum_half += angle.get_data(i,0);
                sum_counts += 1;
                i++;
            }
        j1++;
        }else{
            while(time.get_data(i,0) < env_bot.get_data(j2, 0))
                i++;
            ///суммируем углы между точками огибающей
            while(time.get_data(i,0) < env_top.get_data(j1, 0))
            {//cout << "\t" << time.get_data(i,0) << "\t" << j << "\t" << i;
                sum_half += angle.get_data(i,0);
                sum_counts += 1;
                i++;
            }
        j2++;
        res_str.set_data(0, 0, (env_top.get_data(j1, 0)+env_bot.get_data(j2, 0))/2.0);
        }
        if(sum_counts != 0.0)
            res_str.set_data(0, 1, sum_half/sum_counts);///avg_all_points
        sum_half = 0.0;sum_counts = 0;
        if(res_str.get_data(0, 0) != 0)
            result = result.merge_height(res_str);
    }
    middle_angle = result;
}

matrix oscillation::middle_angle_calc(matrix env_top, matrix env_bot)
{
    matrix result;
    matrix res_str(1,3);///time avg_all_points avg_envelop
    double env_mid_tb = 0.0, env_mid_bt = 0.0, angle_avg = 0.0;
    int i = 0, j = 0, period_counts = 0;

    while( j < env_bot.get_height()-1 &&  j < env_top.get_height()-1)
    {//cout << j << "\t" << env_bot.get_height()-1 << "\t" << env_top.get_data(j, 0) << "\t" << env_bot.get_data(j, 0);
        ///проход по точкам огибающих

        ///вычисляем разницу огибающих
        if( env_top.get_data(j, 0) < env_bot.get_data(j, 0))///верхний раньше нижнего
        {
        ///начинается с верхней огибающей
            //if( j< env_top.get_height()-1)
                env_mid_tb = abs(env_top.get_data(j,1)) - abs(env_bot.get_data(j,1));
            //if( j< env_top.get_height()-1)
                env_mid_bt = -1.0*(abs(env_bot.get_data(j,1)) - abs(env_top.get_data(j+1,1)));
        ///средний угол по точкам
            while(time.get_data(i,0) < env_top.get_data(j, 0))
                i++;
            while(time.get_data(i,0) < env_top.get_data(j+1, 0))
            {//cout << "\t" << time.get_data(i,0) << "\t" << j << "\t" << i;
                angle_avg += angle.get_data(i,0);
                period_counts += 1;
                i++;
            }
        //cout << "\ttop<bot\t" << env_top.get_data(j, 0) << "\t" << angle_avg/period_counts << "\t" << env_mid_tb + env_mid_bt << endl;
        res_str.set_data(0, 0, env_top.get_data(j, 0));///time
        res_str.set_data(0, 1, angle_avg/period_counts);///avg_all_points
        res_str.set_data(0, 2, env_mid_tb + env_mid_bt);///avg_envelop
        period_counts = 0;
        angle_avg = 0.0;
        }else{
        ///начинается с нижней огибающей
            if( j< env_bot.get_height()-1)
                env_mid_bt = -1.0*(abs(env_bot.get_data(j,1)) - abs(env_top.get_data(j,1)));
            if( j< env_top.get_height()-1)
                env_mid_tb = abs(env_top.get_data(j,1)) - abs(env_bot.get_data(j+1,1));
        ///средний угол по точкам
            while(time.get_data(i,0) < env_bot.get_data(j, 0))
                i++;
            while(time.get_data(i,0) < env_bot.get_data(j+1, 0))
            {//cout << "\t" << time.get_data(i,0) << "\t" << j << "\t" << i;
                angle_avg += angle.get_data(i,0);
                period_counts += 1;
                i++;
            }
        //cout << "\tbot>top\t" << env_top.get_data(j, 0) << "\t" << angle_avg/period_counts << "\t" << env_mid_tb + env_mid_bt << endl;
        res_str.set_data(0, 0, env_top.get_data(j, 0));///time
        res_str.set_data(0, 1, angle_avg/period_counts);///avg_all_points
        res_str.set_data(0, 2, env_mid_tb + env_mid_bt);///avg_envelop
        period_counts = 0;
        angle_avg = 0.0;
        }
        j++;
    result = result.merge_height(res_str, result.get_height());
    }
    result.info();
    return result;
};


ostream& operator << (ostream &C, const oscillation &A)
{
    for (int i = 0; i < A.angle.get_height(); i++)
    {
        C << i+1 << "\t\t";
        //C << setw(5) << A.data_source << "\t\t";///убрано для возможности простого чтения файла
        if(A.flow.get_width()==3)
        {
            for (int j = 0; j < A.flow.get_width(); j++)
                C << setw(9) << A.flow.get_data(i,j) << "\t";
        }else{        ///1 Re, 2 M, 6 q
            if(A.flow.get_width()==8)
            {C << setw(9) << A.flow.get_data(i,1) << "\t"
                    << setw(9) << A.flow.get_data(i,2) << "\t"
                        << setw(9) << A.flow.get_data(i,6) << "\t";}
        }
        for (int j = 0; j < A.angle.get_width(); j++)
            C << setw(9) << A.angle.get_data(i,j) << "\t\t";
        for (int j = 0; j < A.dangle.get_width(); j++)
            C << setw(9) << A.dangle.get_data(i,j) << "\t\t";
        C << setw(9) << A.freq.get_data(i,0) << "\t\t";
        C << setw(9) << A.time.get_data(i,0);

        if (i !=  A.angle.get_height()-1)///чтобы не было пустой строки в конце
            C << endl;
    }
    return C;
};

/**oscillation& operator=(const oscillation& A)
{

};*/


///Data Load //Save
void oscillation::Load_Oscillation_Raw_Protocol(const string Data_Way, const string base_file_name)
{cout << "oscillation::Load_Oscillation_Raw_Protocol(const string Data_Way, const string base_file_name)\n";
    matrix m; m.Load_matrix(Data_Way + base_file_name + ".txt");//m.info();//cout << m.get_height() << "\t" << m.get_width() << "\t" << m.get_data(0, 1) << "\t" << m.get_data(m.get_height()-3, 1) << "\t" << m.get_data(m.get_height()-2, 1) << "\t" << m.get_data(m.get_height()-1, 1) << "\n";//cout << m << endl;
        time = m.get_column(0);
        angle = m.get_column(1);
        dangle = this->count_freq(time, angle);
        freq = this->w_harm1();
            freq.erase_column(0);
        matrix zero (m.get_height(), 3); flow = zero;
};

double oscillation::Load_get_Starting_Impolse_Time(string Data_Way, string base_file_name)
{cout << "oscillation::Load_get_Starting_Impolse_Time(string Data_Way, string base_file_name)" << endl;
    ///парсим, ищем в нужном столбце (double)3, берем соответствующее время
    string key = "3";
    string s_buff = "";
    double d_buff = 0;
        ifstream flow_fin(Data_Way + base_file_name + ".ptl");
        while(s_buff != "(сек)")//"N            P0a        синхр           TX          P0b         PC_M           AP           AI         P0_M     t1 (сек)     t2 (сек)"
            flow_fin >> s_buff;
                flow_fin >> s_buff; flow_fin >> s_buff;
        bool flag = true;
        while(flag)
        {
            flow_fin >> s_buff;//cout << s_buff << endl;
            if(s_buff != "______________________________________________________________________________________")
            {
                flow_fin >> s_buff >> s_buff >> d_buff;
                if(d_buff == 3)
                { flag = false; break;}
                for(int i = 0; i < 6; i++)
                    flow_fin >> s_buff;
            }
            else{
                flag = false;
            }
        }
        ///дочитать до искомого элемента времени
        for(int i = 0; i < 6; i++)
                    flow_fin >> s_buff >> d_buff;
        flow_fin.close();
    return d_buff;
}

matrix oscillation::Flow_from_PTL(string Data_Way, string base_file_name)
{cout << "oscillation::Flow_from_PTL(string Data_Way, string base_file_name)\n";
    ifstream fin1(Data_Way + base_file_name + ".ptl");
        string buff = "";
        int b_count = 0;
        int data_height = -2;///????????

    while (b_count < 2)
    {
            fin1 >> buff;
            if(buff == "(сек)")
                b_count++;
    }b_count = 0;
    while(buff!="______________________________________________________________________________________")
    {
        std::getline(fin1,buff);
        //cout << buff << endl;
        data_height++;
    }cout << "data_height = " << data_height << endl;
    fin1.close();

    ifstream fin(Data_Way + base_file_name + ".ptl");

        while (buff != " N             Re            Q     t1 (сек)     t2 (сек)")
            {std::getline(fin,buff);}
    //cout << buff << endl;

    matrix Flow_in(data_height,5);
    double data = 0;
    for(int i = 0; i < data_height; i++)
        for(int j = 0; j < 5; j++)
        {
            fin >> data;
            Flow_in.set_data(i,j, data);
        }
    fin.close();
    //cout << Flow_in << endl;
    for(int i = 0; i < Flow_in.get_height(); i++)
    {
        Flow_in.set_data(i,1,Flow_in.get_data(i,1));///M
        Flow_in.set_data(i,2,Flow_in.get_data(i,2));///Re_l
        Flow_in.set_data(i,3,Flow_in.get_data(i,3)*9.80665);///Q
    }
    return Flow_in;
}

int oscillation::Load_get_Starting_Impolse_Number(string Data_Way, string base_file_name)
{cout << "oscillation::Load_get_Starting_Impolse_Number(string Data_Way, string base_file_name)" << endl;
    ///парсим, ищем в нужном столбце (double)3, берем соответствующее время
    string key = "3";
    string s_buff = "";
    double d_buff = 0;
    int n = 0;
        ifstream flow_fin(Data_Way + base_file_name + ".ptl");cout << base_file_name + ".ptl opened\n";
        while(s_buff != "(сек)")//"N            P0a        синхр           TX          P0b         PC_M           AP           AI         P0_M     t1 (сек)     t2 (сек)"
            flow_fin >> s_buff;
                flow_fin >> s_buff; flow_fin >> s_buff; //cout << s_buff << endl;
        //cout << "**************************************************\n";
        bool flag = true;
        while(flag)
        {
            flow_fin >> s_buff;//cout << s_buff << endl; ///номер
            if(s_buff != "______________________________________________________________________________________")
            {
                flow_fin /*>> s_buff*/ >> d_buff >> s_buff;cout << "S IMP = " << d_buff << endl; ///28.09.19 переставлены столбцы
                n++;
                if(d_buff == 3)
                { flag = false; break;}
                for(int i = 0; i < 8; i++)
                    flow_fin >> s_buff;//cout << s_buff << endl;
            }
            else{
                flag = false;
            }
        }
        flow_fin.close();
    return n;
}

void oscillation::Load_Oscillation_flow_Protocol(string Data_Way, string base_file_name, double t_shift)
{cout << "oscillation::Load_Oscillation_flow_Protocol(string Data_Way, string base_file_name, double t_shift)\n";

    ///поток с дублированием данных
    matrix n_flow(angle.get_height(), 3);
    bool flag = true;

if(FileExists(Data_Way + base_file_name + "_flow.txt")){
    flag = false;
    ///сырые данные потока из чистой матрицы
    matrix buff; buff.Load_matrix(Data_Way + base_file_name + "_flow.txt");//cout << this->Load_get_Starting_Impolse_Number(Data_Way, base_file_name) << endl;

    ///индекс срабатывания синхроимпульса
    int from_index = this->Load_get_Starting_Impolse_Number(Data_Way, base_file_name);///синхроимпульс
        //cout << "flow from index = " << from_index << endl;

    ///выделить указанный промежуток и распространить его на соответсвующие отсчеты угла
    int j = from_index;
    ///double t_shift = angle.get_data(36300, 0) - buff.get_data(j, 4);///по углу реально не с нуля ТУТ ЧУШЬ ///время выхода арретира найти и передавать
    for(int i = 0; i < angle.get_height(); )
    {
        if(angle.get_data(i, 0) < (buff.get_data(j, 4) - buff.get_data(from_index, 4) + t_shift))///записывается один отсчет потока
            while(angle.get_data(i, 0) < (buff.get_data(j, 4) - buff.get_data(from_index, 4) + t_shift ))
            {
                n_flow.set_data(i, 0, buff.get_data(j, 1)); //cout << buff.get_data(j, 1) << "\t";
                    n_flow.set_data(i, 1, buff.get_data(j, 2)); //cout << buff.get_data(j, 2) << "\t";
                        n_flow.set_data(i, 2, buff.get_data(j, 3)); //cout << buff.get_data(j, 3) << "\n";
                i++;
               if(i == angle.get_height()-1)///падало без -1
               { cout << "while break i == angle.get_height();\n"; break;}
            }
        if(n_flow.get_data(i, 1) == 0)
            cout << "angle_index = " <<  i <<"\t" <<  "flow_index = " << j << "\t*******\n";
        if(i == angle.get_height())
                { break;}
        j++;
        if(j == buff.get_height())
        { cout << "while break j == buff.get_height();\n";break;}
    }
    flow = n_flow; n_flow.info();
}
    if(flag)
    {
        matrix F;F.info();
            F = Flow_from_PTL(Data_Way, base_file_name);
        F.Save_matrix(Data_Way + base_file_name + "_flow.txt");
        this->Load_Oscillation_flow_Protocol(Data_Way, base_file_name, t_shift);
    }

};

void oscillation::Load_Oscillation_caliber_angle(string Data_Way, string base_file_name)
{cout << "oscillation::Load_Oscillation_caliber_angle(string Data_Way, string base_file_name)\t";
    string s_buff = "";
    double c_angle = 0;
    ifstream fin(Data_Way + "callibr_angle.txt");
        while(!fin.eof())/// в поисках ключегого слово по всему файлу /// возможно нахождение более одного раза - будет выполнено столько раз
        {
            fin >> s_buff;
                if(s_buff == base_file_name )///file_name -> ключевое слово, за ни следует double сдвижка по времени
                {
                    fin >> c_angle;
                    this->calibr_angle = c_angle;
                }
        }
    fin.close();
    cout << c_angle << endl;
}

void oscillation::Load_Oscillation_Protocol(string Data_Way, string base_file_name)
{
    matrix buff;
        buff.Load_matrix(Data_Way + base_file_name + ".txt");
    ///base
    time = buff.get_column(0);
    angle = buff;
    ///extra
    dangle = this->count_freq(buff.get_column(0), buff.get_column(1));
    freq = this->w_harm1();
        freq.erase_column(0);//freq = freq.merge_width(time, 0);
        angle.erase_column(0);
    ///ifstream fin(Data_Way + base_file_name + "_calibr_angle.txt"); fin >> calibr_angle; fin.close();
       /// this->callibrate_to_radian();
       /// this->callibrate_to_zero();

    flow.Load_matrix(Data_Way + base_file_name + "_flow.txt");
    if(flow.get_height() == 0)
    {
        matrix n_flow(time.get_height(), 3);
        flow = n_flow;
    }
    data_source = base_file_name;
};

void oscillation::Load_Oscillation_Prop_File(string Data_Way, string base_file_name, string mode)
{cout << "oscillation::Load_Oscillation_Prop_File(string Data_Way, string base_file_name)\n";
    matrix B;
        B.Load_matrix(Data_Way + base_file_name + ".txt");
            B.info();
    time = B.get_column(B.get_width()-1);
    freq = B.get_column(B.get_width()-2);
    angle = B.get_column(B.get_width()-4);
    dangle = B.get_column(B.get_width()-3);
    flow = B.get_matrix_part(0,B.get_height(),1,4);

    matrix b_a(B.get_height(),1);
        balance_angle = b_a;
    if(mode == "flow_summary")
    {
       //ifstream fin(Data_Way + base_file_name + );

        //fin.close();
    }
}

void oscillation::Load_Oscillation_Prop_File(string Data_Way, string base_file_name)
{this->Load_Oscillation_Prop_File(Data_Way, base_file_name, "");};

void oscillation::Load_flow_summary(string DataWay, string file_name)
{cout << "Load_flow_summary(string DataWay, string file_name)\n";
    ifstream fin(DataWay + file_name + ".txt");
    string s_buff = "OK";
    double d_buff = 0.0;
    int i_buff = 0;
if(flow.get_data(0,0) != 0 /*&& flow.get_data(flow.get_height(),flow.get_width()) != 0*/)
{
    s_buff = "notOK";
    while( !fin.eof())
        {
            fin >> s_buff; /// номер протокола
                //cout << s_buff << "\t" << s_buff << "\t" << this->get_data_source() << endl;
            if(s_buff == this->get_data_source())
            {
                fin >> s_buff;/// название модели
                    if(this->get_model() == "")
                        this->set_model(s_buff);
                    //cout << s_buff << "\t";
                fin >> i_buff;/// номер первого отсчета потока
                fin >> d_buff;/// число Рейнольдса
                fin >> d_buff;/// число Маха
                fin >> d_buff;/// Тх
                fin >> d_buff;/// скорость потока
                fin >> d_buff;/// скоростной напор в кгс
                fin >> d_buff;/// скоростной напор в Па
                fin >> d_buff;/// плотность газа
            break;
            }
        }
        cout << "WARNING same flow_data\n";
}
    if(s_buff == "OK")
    {
        matrix new_flow(this->flow.get_height(), 8);
            this->set_flow(new_flow);

        while( !fin.eof())
        {
            fin >> s_buff; /// номер протокола
                //cout << s_buff << "\t" << s_buff << "\t" << this->get_data_source() << endl;
            if(s_buff == this->get_data_source())
            {
                fin >> s_buff;/// название модели
                    if(this->get_model() == "")
                        this->set_model(s_buff);
                    cout << s_buff << "\t";
                fin >> i_buff;/// номер первого отсчета потока
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,0,i_buff);
                        cout << i_buff << "\t";
                fin >> d_buff;/// число Рейнольдса
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,1,d_buff);
                        cout << d_buff << "\t";
                fin >> d_buff;/// число Маха
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,2,d_buff);
                        cout << d_buff << "\t";
                fin >> d_buff;/// Тх
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,3,d_buff);
                        cout << d_buff << "\t";
                fin >> d_buff;/// скорость потока
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,4,d_buff);
                        cout << d_buff << "\t";
                fin >> d_buff;/// скоростной напор в кгс
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,5,d_buff);
                        cout << d_buff << "\t";
                fin >> d_buff;/// скоростной напор в Па
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,6,d_buff);
                        cout << d_buff << "\t";
                fin >> d_buff;/// плотность газа
                    for(int i = 0; i < this->get_angle().get_height(); i++)
                        this->set_flow_data(i,7,d_buff);
                        cout << d_buff << "\n";
            break;
            }
        }
    }
    fin.close();
}

void oscillation::Load_model_info(string DataWay, string file_name)
{cout << "Load_model_info(string DataWay, string file_name)\n";
    ifstream fin(DataWay + file_name + ".txt");
    string s_buff = "";
    double d_buff1 = 0.0, d_buff2 = 0.0, d_buff3 = 0.0;
    double key = true; ///продолжать ли поиск?
    while(!fin.eof()&& key)
    {
        fin >> s_buff >> d_buff1 >> d_buff2 >> d_buff3;
        if(s_buff == model_name)
        {
            this->set_model_I(d_buff1);
            this->set_model_S(d_buff2);
            this->set_model_l(d_buff3);
            key = false;
        }
    }
    fin.close();
};

void oscillation::Load_balance_angle_summary(string DataWay, string file_name)
{cout << "Load_middle_angle_summary\t";
if(FileExists(DataWay + file_name + ".txt"))
{
    ifstream fin(DataWay + file_name + ".txt");
    string s_buff = "";
    bool key = false;///использовались ли данные?
    double d_buff1 = 0.0, d_buff2 = 0.0;
    int i = 0;
    while(!fin.eof())
    {
        fin >> s_buff >> d_buff1 >> d_buff2;//cout << s_buff << "\t" << d_buff1 << "\t" << d_buff2 << endl;
        if(s_buff == data_source)
        {
            while( time.get_data(i,0) <= d_buff1 && i < time.get_height()-1)
            {//cout << i << endl;
                balance_angle.set_data(i,0,d_buff2);
                i++;
            }
            key = true;
        }
    }
    fin.close();
    if(!key)///данные не использовались (не было имени data_source)
        {cout << "WARNING\t file balance_angle_summary do not contain " + data_source + " protocol !\n";
    }else{cout << "success\n";}
}else{
cout << "WARNING\t balance_angle_summary file do not exists!\n";
}
}

void oscillation::Save_Oscillation(string Data_Way, string base_file_name)
{
   /* matrix R = time;

    R.merge_width(angle, R.get_width());
    R.merge_width(dangle, R.get_width());
    R.merge_width(flow, R.get_width());

    R.Save_matrix(Data_Way + base_file_name + ".txt");*/
    ofstream fout(Data_Way + base_file_name + ".txt");
    fout << *this;
    fout.close();
}

///
void oscillation::info()
{
    cout << "oscillation object\n";
    cout << "data_source = " << data_source << endl;
    cout << "model_name = " << model_name << "\tI = " << I << "\tS = " << S << "\tl = " << l << endl;
    //cout << "time\t";
      //  time.info();
    cout << "angle\t";
        angle.info();
    cout << "dangle\t";
        dangle.info();
    cout << "flow\t";
        flow.info();
    cout << "calibr_angle = " << calibr_angle << endl;
};



