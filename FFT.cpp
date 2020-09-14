#ifndef FFT_CPP_INCLUDED
#define FFT_CPP_INCLUDED

#include <iostream>
#include <math.h>
#include <fstream>

#include "matrix.h"

using namespace std;

const double TwoPi = 6.283185307179586;
const double Pi = 3.1415926535897932384626433832795;

const double sens_num = 6;
const int time_index = 0;
const int sens_index = 1;

///«а подробност€ми -- “.  ормен, „. Ћейзерсон, –. –ивест,  . Ўтайн, "јлгоритмы. ѕостроение и анализ", ¬торое издание. 2012 г., с. 926-942.
///частота найквиста оцениваетс€ до ~dt/2^N
/// AVal - массив анализируемых данных, Nvl - длина массива должна быть кратна степени 2.
/// FTvl - массив полученных значений, Nft - длина массива должна быть равна Nvl.
void FFTAnalysis(double *AVal, double *FTvl, int Nvl, int Nft)
{
  int i, j, n, m, Mmax, Istp;
  double Tmpr, Tmpi, Wtmp, Theta;
  double Wpr, Wpi, Wr, Wi;
  double *Tmvl;

  n = Nvl * 2; Tmvl = new double[n];

  for (i = 0; i < n; i+=2) {
   Tmvl[i] = 0;
   Tmvl[i+1] = AVal[i/2];
  }

  i = 1; j = 1;
  while (i < n) {
    if (j > i) {
      Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
      Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
    }
    i = i + 2; m = Nvl;
    while ((m >= 2) && (j > m)) {
      j = j - m; m = m >> 1;
    }
    j = j + m;
  }

  Mmax = 2;
  while (n > Mmax) {
    Theta = -TwoPi / Mmax; Wpi = sin(Theta);
    Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
    Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

    while (m < Mmax) {
      i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
      Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
      Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

      while (i < n) {
        j = i + Mmax;
        Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
        Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

        Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
        Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
        i = i + Istp;
      }
    }

    Mmax = Istp;
  }

  for (i = 0; i < Nft; i++) {
    j = i * 2; FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;
  }

  delete []Tmvl;
}


matrix Gauss_Window(matrix F, double sigma)
{
    matrix R(F.get_height(), F.get_width());
    double A = (F.get_height() - 1.0)/2.0, w = 0;

    for(int i = 0; i < R.get_height(); i++)
    {
        w = exp(-pow((i-A)/sigma/A, 2)/2.0);
            R.set_data(i, 0, F.get_data(i,0));
                R.set_data(i, 1, w*F.get_data(i,1));
    }
    return R;
};


double Hann_1_Window(double n, double N)
{
    return(0.5 - 0.5*cos(TwoPi*n/(N-1)));
}

double (*pHann_1_Window)(double, double) = &Hann_1_Window;

double Hamming_Window(double n, double N)
{
    return (0.54 - 0.46*cos(TwoPi*n/(N-1)));
}

double (*pHamming_Window)(double, double) = &Hamming_Window;

void Zero_Addition(double* arr1, double* arr2, int& length, int add)
{cout << "Zero_Addition(double* arr1, double* arr2, int& length, int add)\t" << length << "\tresult length\t";
    if(add > 0)
    { cout << "adding\t";
        double* n_arr1 = new double [length + add];
        double* n_arr2 = new double [length + add];
        for(int i = 0; i < length; i++)
        {
            n_arr1 [i] = arr1[i];
            n_arr2 [i] = arr2[i];
        }
        for(int i = length; i < length + add; i++)
        {
            n_arr1[i] = 0.0;
            n_arr2[i] = 0.0;
        }
        delete [] arr1;
        delete [] arr2;
        arr1 = n_arr1;
        arr2 = n_arr2;

    length += add;
    cout << length << endl;
    }else{cout << "not added\n";}

};

int Spectr_Zero_Addition(int length, int zero_index)
{
    int n_length = 2;
    ///соблюдение длины кратной 2 ///первое увеличение размера
    while(n_length < length)
        n_length *= 2;
    ///размер дополнени€ окна нул€ми ///второе увеличение размера
    if(zero_index > 1)
        n_length *= pow(2, zero_index);
    return n_length;
}

void Get_Spectr(double *f, double* w, double window, int& length, int zero_index)
{cout << "Get_Spectr(double *f, double* w, double window, int& length, int zero_index)\n";
    int n_length = Spectr_Zero_Addition(length, zero_index);
    ///поправка мощности
        double power_k; power_k = n_length/length;
    ///дополнение окна нул€ми
    Zero_Addition(f, w, length, n_length - length);
    FFTAnalysis(f, w, length, length);
    cout << "FFT finished\n";
};

void Get_Spectr(matrix F, matrix& W, double window, int zero_index)
{
    int n_length = Spectr_Zero_Addition(F.get_height(), zero_index);
    ///поправка мощности
        double power_k = n_length/F.get_height();
    matrix n_F(n_length, F.get_width());
    matrix n_W(n_length, F.get_width());
        double w [n_F.get_height()];
        double f [n_F.get_height()];

    for(int j = 0; j < n_F.get_width(); j++)
    {
        ///чтение
        for(int i = 0; i < F.get_height(); i++)
            f[i] = F.get_data(i, j)*power_k;
        for(int i = F.get_height(); i < n_length; i++)
            f[i] = 0;
        ///счет
        FFTAnalysis(f, w, n_length, n_length);
        ///запись
        for(int i = 0; i < n_length; i++)
            n_W.set_data(i, j, w[i]);
    }
    W = n_W;
};

matrix Get_Spectr_column(matrix F, int t_column_index, int f_column_index, int zero_index) ///result freq/power
{//cout << "Get_Spectr_column()\n";//F.info();
    ///подготовка данных
    ///установка длины обработки данных
    int n_length = Spectr_Zero_Addition(F.get_height(), zero_index);
    ///поправка мощности
        double power_k = n_length/F.get_height();

    double* w = new double [n_length];
    double* f = new double [n_length];
        for(int i = 0; i < F.get_height(); i++)
            f[i] = F.get_data(i, f_column_index)*power_k;
        ///заполнение нул€ми
        for(int i = F.get_height(); i < n_length; i++)
            f[i] = 0;
        ///счет
        FFTAnalysis(f, w, n_length, n_length);

    ///результаты
    matrix n_W(n_length, 2);
    ///f=k/T; f=dw*k/N
    //double dw = 1/(F.get_data(0,F.get_height()-1) + (F.get_data(1, f_column_index) - F.get_data(0, f_column_index)));
    double dw = 1.0/fabs(F.get_data(1, t_column_index) - F.get_data(0, t_column_index))/n_length;
    //cout << F.get_data(1, 0) << "\t" <<  F.get_data(0,0) << "\tdw = " << dw << "\tn_length = " << n_length << endl;
        for(int i = 0; i < n_length; i++)
        {
            n_W.set_data(i, 0, dw*i);
            n_W.set_data(i, 1, w[i]);
        }
    delete [] w;
    delete [] f;
    return n_W;
};

matrix Get_Spectr_column_gap(matrix F, int t_column_index, int f_column_index, int from_index, int to_index, int zero_index)
{//cout << "Get_Spectr_column_gap()" << endl;
    matrix f1;
        f1 = F.get_column(t_column_index);
        f1 = f1.merge_width(F.get_column(f_column_index), 1);
        f1 = f1.get_matrix_part(from_index, to_index, 0, 2);
        //cout << f1 << endl;
    return Get_Spectr_column(f1, 0, 1, zero_index);
};

matrix Get_Spectr_column_gap(matrix t, matrix f, int from_index, int to_index, int zero_index)
{//cout << "Get_Spectr_column_gap(matrix t, matrix f, int from_index, int to_index, int zero_index)" << endl;
    matrix f1;
        f1 = t;
        f1 = f1.merge_width(f, f1.get_width());
            f1 = f1.get_matrix_part(from_index, to_index, 0, 2);
    return Get_Spectr_column( f1, 0, 1, zero_index);
};

matrix Get_Spectr_column_gap(matrix F, int t_column_index, int f_column_index, int from_index, int to_index, int zero_index, double sigma)
{
    matrix f1;
        f1 = F.get_column(t_column_index);
        f1 = f1.merge_width(F.get_column(f_column_index), 1);
        f1 = f1.get_matrix_part(from_index, to_index, 0, 2);
    return Get_Spectr_column( Gauss_Window(f1, sigma), 0, 1, zero_index);
}

matrix Get_Spectr_column_gap(matrix t, matrix f, int from_index, int to_index, int zero_index, double sigma)
{
    matrix f1;
        f1 = t;
        f1 = f1.merge_width(f, f1.get_width());
            f1 = f1.get_matrix_part(from_index, to_index, 0, 2);
            f1 = Gauss_Window(f1, sigma);
    return Get_Spectr_column( f1, 0, 1, zero_index);
};

matrix Get_Spectr_max_w(matrix S)
{//cout << "Get_Spectr_max_w(matrix S)\n";
    matrix R(1,2);
    double m = 0; int index_m = 0;
    for(int i = 0; i < (int)S.get_height()/2; i++)///проход по частотам
        if(S.get_data(i,1) > m)///если мощность больше заданной
        {
            m = S.get_data(i,1);///то она максимальна€, зададим еЄ
            index_m = i; ///запишем где нашли
        }
    R.set_data(0,1, S.get_data(index_m,0));/// значение частоты будет во второй €чейке
    R.set_data(0,0, S.get_data(index_m,1));/// значение мощнсти в первой
    return R;
};


matrix Front_Get_Main_freq(matrix Signal, matrix env_top, matrix env_bot)
{//cout << "Front_Get_Main_freq(matrix Signal, matrix env_top, matrix env_bot)\n";
matrix R(1,2);
    double T0 = 0;
    int sum_counter = 0;

    double t_min, t_max;
        t_min = Signal.get_data(0,0);
        t_max = Signal.get_data(Signal.get_height()-1,0);
    int k = 0;
    for(int i  = 1; i < env_top.get_height(); i++)
        if( env_top.get_data(i-1,0) >= t_min && env_top.get_data(i,0) <= t_max )
        {//cout << i<< endl;
            T0 += (env_top.get_data(i,0) - env_top.get_data(i-1,0));
            k++;
        }
    T0 = T0/k;
    cout << t_min <<"\t" << t_max << "\t" << T0 << "\t" << 1.0/T0 << endl;
    R.set_data(0,0,0);
    R.set_data(0,1, 1.0/T0);
return R;
}

matrix Get_Freq_Window_FFT(matrix D, int step_size, int window_size, int zero_index)
{
            matrix f_m, spectr, f;
            //int window_counter = 0;
            D.info();
            //int k;
                //cin >> k;
            for(int index_h = 0; index_h < D.get_height()-window_size;)
            {cout << "FFT index_h = " << index_h << "\t" << (double)index_h/D.get_height()*100 << "%\t";
                spectr = Get_Spectr_column_gap(D.get_column(0), D.get_column(1), index_h, index_h + window_size, zero_index);
                    f = Get_Spectr_max_w(spectr);///плотность и наиболее плотна€ частота между index_h и index_h + window_size
                    f.set_data(0,0, D.get_column(0).get_data(index_h + (int)window_size/2,0));///врем€ в середине участка и частота главной гармоники ///убрано!!!!!!!
                f_m = f_m.merge_height(f, f_m.get_height());///собираем полученные частоты
                index_h +=step_size;
                cout /*<< f.get_data(0,0)*/ << "\tfinished\n";
            }
            return f_m;
}

#endif // FFT_CPP_INCLUDED
