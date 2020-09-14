#include <iostream>

#include <cmath>

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <string>

#include "matrix.h"
#include "Approximation.h"

#include "Rot_Mx.h"

using namespace std;

const double PI = 3.14159265359;

Rot_Mx::Rot_Mx()
{
    matrix t(0,1);
    //time = t;
    matrix v(0,1);
    voltage = v;
    matrix f(0, 1);
    freq = f;
    time_shift = 0;
    holl_conf = "";

    matrix veloc(0,1);
    velocity = veloc;
};

/*void Rot_Mx::Calc_velocity()
{
    velocity = this->Get_Calc_velocity();
};*/

void Rot_Mx::set_holl_conf(string h_c)
    {holl_conf = h_c;};

matrix Rot_Mx::Get_Calc_velocity()
{cout << "matrix Rot_Mx::Get_Calc_velocity()\n";
    double LIM;
    ///отчеты ацп, положительное напряжение = отрец отчет
    LIM = -200; ///изменил знаки сравнения

    bool flag = 1; ///не находили
    double time0 = 0, time1 = 0, time_c = 0, time_c_prev = 0; ///первый, второй, между ними, и между ними предыдущий
    int i_prev = 0;///номер отсчета прошлого пролета магнита
    int prev_begin_num = 0;
    //double dt1 = 0, dt2 = 0;

    double* V_arr = new double[voltage.get_height()];
    double* t_arr = new double [voltage.get_height()];

    int i = 0, j = 0;
    while( i < voltage.get_height()-2)
    {///cout << "while\t" << i << endl;

        while( (voltage.get_data(i,1) >= LIM) && (i < voltage.get_height()-1) )
            i++;
        ///нашли пик и записали
            time0 = voltage.get_data(i,0);
        while((voltage.get_data(i,1) < LIM) && (i < voltage.get_height()-1))
            i++;
        ///нашли конец пика и записали
            time1 = voltage.get_data(i-1,0);///сдвинул на 1
        //if(time0 < 45)cout << "time0 =" << time0 << "\ttime1 =" << time1 << "\t" << i << endl;
        ///считаем, что истинное время пролета магнита посередине
            time_c = (time1 + time0)/2;
        //cout << "timec =" << time_c << "\ttime_c_prev\t" << time_c_prev << endl;
        if(/*time_c_prev != 0*/ i>0)
        {
            V_arr[j] = 1/(time_c - time_c_prev); ///частота 1/(время между пролетами магнита)
            t_arr[j] = time_c - (time_c-time_c_prev)/2.0;//time.get_data(i ,0); ///сдвинули на половину времени пролета
        }
        time_c_prev = time_c;
        i_prev = i;
        j++;
            //for(int j = prev_begin_num; j < i; j++ )  {V.set_data(j, 0, 1/(time_c - time_c_prev)*length*PI/5.0);  } ///rad/sec
    }

    int length = 0;
    for(int k = 1; k < voltage.get_height(); k++)
        if(t_arr[k] != 0)
            length++;
//length = 20000;
///cout << "!!!!!!!!!!!!!\tlength\t" << length << endl;
    ///пересчет по количеству магнитов
    double v_coeff = 1;
    if(holl_conf == "1_0_0_0_0_1_0_0_0_0")
        v_coeff = 0.5;

    matrix R(length, 2);
    for(int i = 0, j = 0; i < voltage.get_height(), j < length; i++, j++)
        if(t_arr[i] != 0)
        {
            R.set_data(j,0, t_arr[i]);
            R.set_data(j,1, V_arr[i]*v_coeff);
        }

    delete [] V_arr;
    delete [] t_arr;
    R.info();

   // R = smooth_matrix_column(R, 3, 1);

   /// R.erase_row(0);
    return R;
};

matrix Rot_Mx::Get_Calc_velocity_slope(int window_range)
{
    //cout << this->velocity << endl;
    matrix Slope(this->velocity.get_height(), 2);///краевые выбросили
    double sl = 0;
    matrix Part(window_range, 2);
    int inter;
    for(int i = 0; i < this->velocity.get_height() - window_range; i++)
    {//cout << this->velocity.get_height() << "\t" << i << endl;
        inter = window_range/2;

        for(int j = 0; j < window_range; j++)
        { //cout << "i + j = " << i + j << "\t";
            Part.set_data(j, 0, this->velocity.get_data(i+j, 0));
            Part.set_data(j, 1, this->velocity.get_data(i+j, 1));
            //cout << this->time.get_data(i+j, 0) << "\t" << this->velocity.get_data(i+j, 0) << endl;
         }
        //cout << Part << endl;
        Find_Incline_Coefficient(Part, &sl, window_range);
        //cout << "\nslope = " << sl  << "\t" << Part.get_data(inter+1,0) << endl;

        /*if((int)window_range/2 == (double)window_range/2)
            {Slope.set_data(i, 0, (Part.get_data((int)window_range/2, 0) + Part.get_data((int)window_range/2 + 1), 0)/2.0);}
            else
            {Slope.set_data(i, 0, Part.get_data((int)window_range/2+1), 0);}*/
        Slope.set_data(i, 0, Part.get_data(inter+1,0));
        Slope.set_data(i, 1, sl);
    }
    return Slope;
}


matrix Rot_Mx::Get_Avg_velocity(Rot_Mx A, double from_v)
{cout << "Rot_Mx::Get_Avg_velocity(Rot_Mx A, double from_v)\t";
    ///сдвинули по времени
    velocity.column_shift(0, get_time_shift(velocity, from_v)); ///i1
    (A.velocity).column_shift(0, get_time_shift(A.velocity, from_v));///i2
    cout << "time shifted\n";
    cout << "************" <<  endl;

    velocity.info();
    (A.velocity).info();


    ///сдвинули счетчики в точку начала осреднения (точка "перегиба")
    int i1 = 0, i2 = 0;
    cout << i1 << "\t" << velocity.get_height() << "\t" << velocity.get_data(i1,1) << "\t" <<  from_v << endl;
        while(velocity.get_data(i1, 1) > from_v)
        {//cout << i1 << "\t" << velocity.get_height() << "\t" << velocity.get_data(i1,1) << "\t" <<  from_v << endl;
            i1++;
        }
        while(A.velocity.get_data(i2,1) > from_v)
        {//cout << i2 << "\t" << A.velocity.get_height() << "\t" << A.velocity.get_data(i2,1) << "\t" <<  from_v << endl;
            i2++;
        }
    cout << i1 << "\t" << i2 << "\tshifted\n";

    matrix R(velocity.get_height() + A.velocity.get_height() - i1 - i2, 2);///обрезали до пересечения
    R.info();
    int R_size = 0, i3 = 0;
        ///уже сдвинуты куда нужно
        R.set_data(0, 0, velocity.get_data(i1, 0));
        R.set_data(0, 1, (velocity.get_data(i1, 1) + A.velocity.get_data(i2, 1))/2.0);
        R_size++;

    i1++, i2++; i3++;
    double i10 = i1, i20 = i2;

    cout << "\nSTAGE 1\n";
    ///1. вписать имена в R ///или можно просто вписать и отсортировать (гораздо дольше но проще записать)
    //for(int i = 0; i < min(A.velocity.get_height(), velocity.get_height()); i++)
    for(i3 = 1; i3 < R.get_height()-1; i3++)
    {//cout << i3 << "\t";
        if(velocity.get_data(i1,0) <= A.velocity.get_data(i2, 0))
        {//cout << "if1\t" << i3 << "\t" << R.get_height()-1 << "\t" << i1 << "\t" << velocity.get_height() << "\t" << i2 << "\t" << A.velocity.get_height() << endl;
            R.set_data(i3, 0, velocity.get_data(i1,0));
            i1++;

        }else{//cout << "if2\t" << i3 << "\t" << R.get_height()-1 << "\t" << i1 << "\t" << velocity.get_height() << "\t" << i2 << "\t" << A.velocity.get_height() << endl;
            R.set_data(i3, 0, A.velocity.get_data(i2,0));
            i2++;
        }
        ///остановка
        if(i1 >= velocity.get_height()-1)
        {
            R_size = i3;
            i3 = R.get_height();
        }
        if(i2 >= A.velocity.get_height()-1 )
        {
            R_size = i3;
            i3 = R.get_height();
        }
    }

    cout << "\nSTAGE 2\n";
    ///2. установили верный размер
    for(int k = R.get_height()-1; k >= 1; k--)
    {//cout << R.get_height() << "\t" << R_size << "\t" << k << endl;
        if(R.get_data(k, 0) == R.get_data(k-1, 0))
            R_size--;
        if(R.get_data(k, 0) != R.get_data(k-1, 0))
            k = -1;
    }
    cout << "R height&size = " << R.get_height() << "\t" << R_size << endl;

    for(int k = R.get_height(); k > R_size; k--)
    {//cout << k << endl;
        R.erase_row(R.get_height()-1);
    }
    R.info();

    cout << "\nSTAGE 3\n";
    ///3.
    i1 = i10; i2 = i20;
    matrix m(2,2);
    double sl = 0, v_old = 0, v_new = 0;
    for(int i = 1; i < /*R.get_height()*/305; i++)
    {
        if(velocity.get_data(i1, 0) == R.get_data(i,0))
        {///значит ищем добавочную точку на A
            cout << "if1\t" <<  i << "\t";
            m.set_data(0, 0, A.velocity.get_data(i2 - 1, 0));
            m.set_data(0, 0, A.velocity.get_data(i2 - 1, 1));
            m.set_data(0, 0, A.velocity.get_data(i2, 0));
            m.set_data(0, 0, A.velocity.get_data(i2, 1));
                Find_Incline_Coefficient(m, &sl, m.get_height());
            //R.velocity.set_data(i, 1, velocity.get_data(i2, 1) + sl*(velocity.get_data(i2, 0) - velocity.get_data(i2 - 1, 0));
            v_old = velocity.get_data(i1, 1);
            v_new = velocity.get_data(i2, 1) + sl*(velocity.get_data(i2, 0) - velocity.get_data(i2 - 1, 0));
            i2++;
        }else{
        ///аналогично на v
            cout << "if2\t" <<  i << "\t";
            m.set_data(0, 0, velocity.get_data(i1 - 1, 0));
            m.set_data(0, 0, velocity.get_data(i1 - 1, 1));
            m.set_data(0, 0, velocity.get_data(i1, 0));
            m.set_data(0, 0, velocity.get_data(i1, 1));
                Find_Incline_Coefficient(m, &sl, m.get_height());
            ///R.velocity.set_data(i, 1, velocity.get_data(i1, 1) + sl*(velocity.get_data(i1, 0) - velocity.get_data(i1 - 1, 0));
            v_old = A.velocity.get_data(i2, 1);
            v_new = velocity.get_data(i1, 1) + sl*(velocity.get_data(i1, 0) - velocity.get_data(i1 - 1, 0));
            i1++;
        }

        R.set_data(i, 1, (v_old + v_new)/2.0);
        cout << R.get_data(i, 1) << "\n";
    }

    return R;
}




void Rot_Mx::Load_Rot_Mx(string file_name)
{
    matrix M;
    M.Load_matrix(file_name);

    //time = M.get_column(0);
    voltage = M.get_column(1);

    matrix m(0, 1);
    velocity = m;
    freq = m;
};

void Rot_Mx::Load_Rot_Mx_velocity(string file_name)
{
    matrix M;
    M.Load_matrix(file_name);

    //time = M.get_column(0);
    velocity = M;

    matrix m(0, 1);
    voltage = m;
    freq = m;
};


void Rot_Mx::Save_Rot_Mx(string file_name)
{
    ///matrix R;
    //R = time;
    ///R = R.merge_width( velocity, R.get_width());
    ///R = R.merge_width( w_slope.get_column(1), R.get_width());
    ///R = R.merge_width( w_2slope.get_column(1), R.get_width());
    ///R.Save_matrix(file_name);
    velocity.Save_matrix(file_name);
};


void Rot_Mx::info()
{
    cout << "Rot_Mx object\n";
    //cout << "time\t";
    //time.info();
    cout << "voltage\t";
    voltage.info();
    cout << "freq\t";
    freq.info();
    cout << "velocity\t";
    velocity.info();
}


double get_time_shift(matrix Data, double f)///f частота сравнения
{///Data -> time/w
    int i = 0;
    while( Data.get_data(i,1) >= f)
        i++;//cout << "Data_Row\t" << Data.get_row(i-1) << Data.get_row(i) << "\n" << Data.get_row(i-1) << endl;
    ///можно точнее (линейно) аппроксимируя
    return (Data.get_data(i,0));
}

matrix Rot_Mx_time(string Data_Way, string data_file_name, double ADC_freq)
{
    matrix ADC_data;
    ///Загрузка данных напряжения с АЦП в отчетах (частота опроса ADC_freq)
        ADC_data.Load_matrix(Data_Way + data_file_name + ".txt");
    ///Построение времени по частоте опроса ADC_freq
    matrix time(ADC_data.get_height(),1);
        for (int i = 0; i < ADC_data.get_height(); i++)
            time.set_data(i, 0, (double)i/ADC_freq);
    ADC_data = ADC_data.merge_width(time, 0);
    return ADC_data;
}

void Rot_Mx_test(string Data_Way, string data_file_name, string h_c, double ADC_freq, int slope_window_range, int slope2_window_range)
{cout << "Rot_Mx_test Started\t";
    matrix ADC_data;
    ///Загрузка данных напряжения с АЦП в отчетах (частота опроса ADC_freq)
        ADC_data.Load_matrix(Data_Way + data_file_name + ".txt");

    ///Построение времени по частоте опроса ADC_freq
    matrix time(ADC_data.get_height(),1);
        for (int i = 0; i < ADC_data.get_height(); i++)
            time.set_data(i, 0, (double)i/ADC_freq);
            ///(t.merge_width(ADC_data, t.get_width())).Save_matrix(Data_Way + data_file_name + "_t.txt");
    ADC_data = ADC_data.merge_width(time, 0);

    ///построение объекта Rot_Mx
    Rot_Mx R;
        R.voltage = ADC_data;
            ADC_data.Save_matrix(Data_Way + data_file_name + "_data.txt");

    R.set_holl_conf(h_c);
    ///w ///velocity Nx2
    R.velocity = R.Get_Calc_velocity();
        ///удаление счетных ошибок
        R.velocity.erase_row(0);
        ///сохранение
        R.velocity.Save_matrix(Data_Way + data_file_name + "_w.txt");

    matrix W;
    W = R.velocity;

    ///dw/dt
    matrix V_slope;
        V_slope = Find_Slope_of_Matrix_Column(R.velocity, 0, 1, slope_window_range);
        (V_slope).Save_matrix(Data_Way + data_file_name + "_dw.txt");
        ///W = W.merge_width(V_slope.get_column(1), W.get_width());

    ///d2w/dt2
    matrix V_slope_2;
        V_slope_2 = Find_Slope_of_Matrix_Column(V_slope, 0, 1, slope2_window_range/*(int)V_slope.get_height()/10*/);
        (V_slope_2).Save_matrix(Data_Way + data_file_name + "_d2w.txt");
        ///W = W.merge_width(V_slope_2.get_column(1), W.get_width());

    ///W.Save_matrix(Data_Way + data_file_name + "_w.txt");
cout << "finished\n";
}
