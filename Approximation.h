#ifndef APPROXIMATION_H_INCLUDED
#define APPROXIMATION_H_INCLUDED

#include "matrix.h"

template<typename T>
void bubble_sort(T array[], std::size_t size);

///сортировка столбцов матрицы
int get_max_element_index_from_matrix_column(matrix m, int column_index, int index_from, int index_to);
double get_max_element_from_matrix_column(matrix m, int column_index, int index_from, int index_to);
double get_max_element_from_matrix_column(matrix m, int column_index);

matrix sort_matrix_column(matrix* );//!!
matrix sort_matrix_column(matrix* D, int c_index);

double get_median(double* arr, int length); /// медиана массива
matrix get_median_matrix_column(matrix* M); /// медиана столбцов матрицы

///HISTOGRAM
matrix Get_Envelop_Sum_Histogram(matrix A, int bars_num);///first column summ

///APPROXIMATION
void Find_Incline_Coefficient(matrix m, double* incline, int n); /// m[][0] = x, m[][1] = y; /// f = a*x+b ///наклон
void Find_Incline_Shift_Coefficient(matrix m, double* incline, double* shift, int n); /// m[][0] = x, m[][1] = y; ///наклон и свдиг по y

matrix Find_Slope_of_Matrix_Column(matrix m, int index_x, int index_y, int window_range); ///y(x) = a*x + b /// return(matrix R(xi, ai)) ///R.height = m.height - window_range

int get_plateau_backward(matrix m, int x_inde, int y_ind, double plateau_condition);///выход функции "на полочку"
int get_plateau_forward(matrix m, int x_inde, int y_ind, double plateau_condition);

///filters
matrix Get_Envelope_both1(matrix A, matrix dA);
matrix Get_Envelope_top(matrix A, matrix dA);
    matrix Get_Envelope_top(matrix t, matrix f, matrix df);
    matrix Get_Envelope_top_hard(matrix A, matrix dA, double level);
matrix Get_Envelope_top(matrix A, int index_y, int index_x, matrix dA);
    matrix Get_Envelope_top_hard(matrix A, int index_y, int index_x, matrix dA, double level);///level - уровень меньше которого отсеиваютс€ возможные кандидаты (например уровень вокруг которого совершаютс€ колебани€)
matrix Get_Envelope_top(matrix y_column, matrix x_column, matrix dy_column);/// f(x)//x//df(x)/dx
    matrix Get_Envelope_bot(matrix y_column, matrix x_column, matrix dy_column);/// f(x)//x//df(x)/dx ///invert comparation
matrix Get_Envelope_top_hard(matrix y_column, matrix x_column, matrix dy_column, double level);

matrix Get_Envelop_frq_top(matrix x, matrix y, double freq, double ambit);///ambit дол€ окна от периода /// коррекци€ поиска огибающей окном в зависимости от частоты и пар-ра
    matrix Get_Envelop_frq_bot(matrix x, matrix y, double freq, double ambit);
matrix Get_Envelope_bot(matrix A, matrix dA);
    matrix Get_Envelope_bot(matrix A, int index_y, int index_x, matrix dA);
matrix Get_Envelope_topbot(int index_govern, matrix env_t, matrix env_b);///верхн€€ и нижн€€ огибающие объеденины ///столбцы данных, один из них определ€ющий, аналогично времени

void Window_Sigma_Filter(double* arr, const int arr_length, int window_range, const double sigma_level);
void Filter_Outlier(double* arr, int from, int arr_length, double in_coeff, double out_coeff); /// фильтрование выбросов через квартили///quartiles
void Filter_Outlier1(double* in_arr, int arr_length, double in_coeff, double out_coeff);
void smooth(double *input, double *output, int n, int window);
matrix smooth_matrix_column(matrix A, int window_length, int column_index);
matrix smooth_matrix_height(matrix A, int window_length);
    matrix smooth_matrix_height(matrix A, matrix window_lengthes); ///special window_length for each column

///interpolation
/*//Ќ≈ «ј¬≈–Ў≈Ќќ***************
void Newton_Polinom(double* arr, int arr_length, double* coeff, int polyno_power);
double Polynom_Lagrange(double argument, double* arr_points, int poly_index, int points_amount);
double Get_Lagrange(double argument, double* arr_func, double* arr_points, int points_amount);
double Get_Simple_MA(double* arr, int arr_length, int window_range);

void Window(double* arr, int arr_length, int window_length);

double Standart_Deviation(double* arr, double* sigma,  int window_length );
double Student_Coeff(double P, int f);///подгружать таблицу?

double part_dispertion(double* arr, int data_length);
double part_ref_dispertion(double* arr, int data_length);

void arr_dispertion(double* arr, double* dispertion, int arr_length, int window_range);
void arr_ref_dispertion(double* arr, double* dispertion, int arr_length, int window_range);
///Ќ≈ «ј¬≈–ЎЌќ*****************/


///DEREVATIVES
//double* derevative_3_dot(double* y, double* x);
matrix derevative_3_dot(matrix Y, int y_index, int x_index);///дл€ равномерной сетки
matrix derevative_3_dot_h(matrix Y, int y_index, int x_index);///дл€ равномерной сетки
matrix derevative_3_dot(matrix y, matrix x);///дл€ равномерной сетки

///MATH FUNC
double factorial2(double arg);///четный или не четный факториал (шаг2)
double Legander_K(double arg, int degree);///значение полинома лежандра от аргумента arg степени degree
void equation_2order(const double a, const double b, const double c, double& root1Re, double& root1Im, double& root2Re, double& root2Im);
double matrix_func(matrix m, double arg, int index_arg, int index_func); ///поточечный возврат функции /// добавить интерпол€цию
int sign(double a);

#endif // APPROXIMATION_H_INCLUDED
