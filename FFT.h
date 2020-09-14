#ifndef FFT_H_INCLUDED
#define FFT_H_INCLUDED

void FFTAnalysis(double *AVal, double *FTvl, int Nvl, int Nft);

matrix Gauss_Window(matrix F, double sigma);///возвращает сигнал умноженный на гаусс

void Zero_Addition(double* arr1, double arr2, int length, int add);

void Get_Spectr(double *f, double* w, double , int& length, int zero_inde);

void Get_Spectr(matrix F, matrix& W, double window, int zero_index);
matrix Get_Spectr_column(matrix F, int t_column_index, int f_column_index, int zero_index);///равномерный шаг по вермени ///result freq/power
    matrix Get_Spectr_column_gap(matrix F, int t_column_index, int f_column_index, int from_index, int to_index, int zero_index);
    matrix Get_Spectr_column_gap(matrix t, matrix f, int from_index, int to_index, int zero_index);
    matrix Get_Spectr_column_gap(matrix F, int t_column_index, int f_column_index, int from_index, int to_index, int zero_index, double sigma);
    matrix Get_Spectr_column_gap(matrix t, matrix f, int from_index, int to_index, int zero_index, double sigma);

matrix Get_Spectr_max_w(matrix S);

matrix Front_Get_Main_freq(matrix Signal, matrix env_top, matrix env_bot);

double Hann_1_Window(double n, double N);

double Hamming_Window(double n, double N);

matrix Get_Freq_Window_FFT(matrix D, int step_size, int window_size, int zero_index);

#endif // FFT_H_INCLUDED
