#include <iostream>
#include <cmath>

#include "Approximation.h"
#include "matrix.h"

#include <fstream>
using namespace std;


#include <cstddef>
#include <utility>


template<typename T>
void bubble_sort(T array[], std::size_t size)
{
    for (std::size_t idx_i = 0; idx_i < size - 1; idx_i++)
    {
        for (std::size_t idx_j = 0; idx_j < size - idx_i - 1; idx_j++)
        {
            if (array[idx_j + 1] < array[idx_j])
            {
                std::swap(array[idx_j], array[idx_j + 1]);
            }
        }
    }
}

int get_max_element_index_from_matrix_column(matrix m, int column_index, int index_from, int index_to)
{
    double max_elem = m.get_data(index_from, column_index);
    int max_elem_index = 0;
    for(int i = index_from; i < index_to; i++)
        if(max_elem < m.get_data(i, column_index))
        {
            max_elem = m.get_data(i, column_index);
            max_elem_index = i;
        }
    return max_elem_index;
};

double get_max_element_from_matrix_column(matrix m, int column_index, int index_from, int index_to)
{
    double max_elem = m.get_data(index_from, column_index);
    for(int i = index_from; i < index_to; i++)
        if(max_elem < m.get_data(i, column_index))
            max_elem = m.get_data(i, column_index);
    return max_elem;
};

double get_max_element_from_matrix_column(matrix m, int column_index)
    { return get_max_element_from_matrix_column(m, column_index, 0 , m.get_height()); }


matrix sort_matrix_column(matrix* D)
{
    matrix R; R = *D;
    if( R.get_height() >= 1 && R.get_width() >= 0)
    {
        double arr[R.get_height()];

        for(int j = 0; j < R.get_width(); j++)
        {
            ///������
            for(int i = 0 ; i < R.get_height(); i++)
                arr[i] = R.get_data(i,j);
            ///����
            bubble_sort(arr, R.get_height());

            ///������
            for(int i = 0 ; i < R.get_height(); i++)
                R.set_data(i, j, arr[i]);
        }
    return R;
    }else{cout << "Wrong sizes sort_matrix_column(matrix* ) ! \n"; return *D;}
}

matrix sort_matrix_column(matrix* D, int c_index)
{
    matrix R; R = *D;
    if( R.get_height() >= 1 && R.get_width() >= 0)
    {//cout << c_index << endl;
        double arr[R.get_height()];
            ///������
            for(int i = 0 ; i < R.get_height(); i++)
                arr[i] = R.get_data(i,c_index);
            ///����
                bubble_sort(arr, R.get_height());
            ///������
            for(int i = 0 ; i < R.get_height(); i++)
                R.set_data(i, c_index, arr[i]);
    return R;
    }else{cout << "Wrong sizes sort_matrix_column(matrix* ) ! \n"; return *D;}
}

double get_median(double* arr, int length) /// ������� �������
{
    double s_arr [length];
        for(int i = 0; i < length; i++)
            s_arr[i] = arr[i];
        bubble_sort(s_arr, length);

    if(length%2 == 0)
        return (double)(s_arr[length/2] + s_arr[length/2 + 1])/2;///���
    return s_arr[length/2 + 1];///�����
}

matrix get_median_matrix_column(matrix* M) /// ������� �������� �������
{
    matrix median(1, M->get_width());

    if(M->get_height()%2 == 0)
    {///���
        matrix s_M = sort_matrix_column(M);
        for(int j = 0; j < M->get_width(); j++)
            median.set_data(0, j, (s_M.get_data(M->get_height()/2, j) + s_M.get_data(M->get_height()/2 + 1, j))/2);
    }else
    {///�����
        matrix s_M = sort_matrix_column(M);
        for(int j = 0; j < M->get_width(); j++)
            median.set_data(0, j, s_M.get_data(M->get_height()/2 + 1, j));
    }
    return median;
}

void Find_Incline_Coefficient(matrix m, double* incline, int n)/// Ap.width = 2
{//cout << "matrix in =\n" << m << endl;
    double sumx = 0;
    double sumy = 0;
    double sumx2 = 0;
    double sumxy = 0;
    for (int i = 0; i<n; i++)
    {//cout << "incline\t" << i << "/" << n << endl;
        sumx += m.get_data(i, 0);
        sumy += m.get_data(i, 1);
        sumx2 += m.get_data(i, 0) * m.get_data(i, 0);
        sumxy += m.get_data(i, 0) * m.get_data(i, 1);
    }
  *incline = (n*sumxy - (sumx*sumy)) / (n*sumx2 - sumx*sumx);
};

void Find_Incline_Shift_Coefficient(matrix m, double* incline, double* shift, int n)
{//cout << "matrix in =\n" << m << endl;
    if(m.get_height() >= n)
    {double sumx = 0;
    double sumy = 0;
    double sumx2 = 0;
    double sumxy = 0;
    for (int i = 0; i<n; i++)
    {
        sumx += m.get_data(i, 0);
        sumy += m.get_data(i, 1);
        sumx2 += m.get_data(i, 0) * m.get_data(i, 0);
        sumxy += m.get_data(i, 0) * m.get_data(i, 1);
    }
  *incline = (n*sumxy - (sumx*sumy)) / (n*sumx2 - sumx*sumx);
  *shift = (sumy - *incline*sumx) / n;
  }else{
  *incline = 0;
  *shift = 0;
  }
}

matrix Find_Slope_of_Matrix_Column(matrix m, int index_x, int index_y, int window_range)
{
    matrix Slope(m.get_height()-window_range, 2);///� ������ �������� ���������� �� ��������������
        double sl = 0;
    matrix Part(window_range, 2);
    int inter;
    for(int i = 0; i < m.get_height() - window_range; i++)
    {//cout << this->velocity.get_height() << "\t" << i << endl;
        inter = window_range/2;
        for(int j = 0; j < window_range; j++)
        { //cout << "i + j = " << i + j << "\t";
            Part.set_data(j, 0, m.get_data(i+j, index_x));
            Part.set_data(j, 1, m.get_data(i+j, index_y));
        }
        Find_Incline_Coefficient(Part, &sl, window_range);
        Slope.set_data(i, 0, Part.get_data(inter,0));
        Slope.set_data(i, 1, sl);
    }
    return Slope;
};

int get_plateau_backward(matrix m, int x_inde, int y_ind, double max_plateau_condition)
{cout << "get_plateau_backward(matrix m, int x_inde, int y_ind, double max_plateau_condition)\n";
    double avg = 0;
    int back_range = 2;
    double *arr_of_avg = new double [m.get_height()];
    /// ����� �  ����� ������ ���� ����� ����������, ���� ��� ��� - ��� �����.
    ///�������� "�����������" ���������� �������� plateau_condition
    /// plateau_condition ��� ������� ����� ������� ���������� �� ������ ���� � ����� ����������� ������, ��������������� �� ������� ������� ���������
    double plateau_condition = 0;
    while(plateau_condition <= max_plateau_condition && (m.get_height()-back_range-1) >=0)
    {
        for (int i = m.get_height()-1; i >= m.get_height()-back_range; i--) ///���� ����� � �����
        {
            avg += m.get_data(i,y_ind);
           // cout << i << "\t" << m.get_data(i,y_ind) << "\t" << avg << endl;
        }
        avg = avg/(back_range);
        arr_of_avg[back_range] = avg;
          //  cout << "AVG = " <<  avg <<"\tback_range = " << back_range << endl;
        plateau_condition = fabs((fabs(avg) - fabs(m.get_data(m.get_height()-back_range-1,y_ind))))/fabs(avg);
        //cout << back_range << "\tavg = " << 1 << "\tnew_dot = " << fabs(m.get_data(m.get_height()-back_range-1,y_ind))/fabs(avg) << "\t" << plateau_condition << "\t" << m.get_height()-back_range << endl;
        back_range++;
        avg = 0;//cin >> avg;
    }
    back_range--;
    ///m.get_height()-back_range - ������ �����
    ///cout << "PLATEAU of " << max_plateau_condition << " level FROM " << m.get_height()-back_range << "\t = " << m.get_data(m.get_height()-back_range,y_ind) << endl;
    delete [] arr_of_avg;
    return back_range;
};

int get_plateau_forward(matrix m, int x_inde, int y_ind, double max_plateau_condition)
{cout << "get_plateau_forward(matrix m, int x_inde, int y_ind, double max_plateau_condition)\n";
    double avg = 0;
    int for_range = 2;
    double *arr_of_avg = new double [m.get_height()];
    /// ����� �  ����� ������ ���� ����� ����������, ���� ��� ��� - ��� �����.
    ///�������� "�����������" ���������� �������� plateau_condition
    /// plateau_condition ��� ������� ����� ������� ���������� �� ������ ���� � ����� ����������� ������, ��������������� �� ������� ������� ���������
    double plateau_condition = 0;
    while(plateau_condition <= max_plateau_condition && (for_range <= m.get_height()-1))
    {
        for (int i = 0; i < for_range; i++) ///���� ����� � �����
        {
            avg += m.get_data(i,y_ind);
         //   cout << i << "\t" << m.get_data(i,y_ind) << "\t" << avg << endl;
        }
        avg = avg/(for_range);
        arr_of_avg[for_range] = avg;
        //    cout << "AVG = " <<  avg <<"\tfor_range = " << for_range << endl;
        plateau_condition = fabs((fabs(avg) - fabs(m.get_data(for_range+1,y_ind))))/fabs(avg);
        //cout << for_range << "\tavg = " << 1 << "\tnew_dot = " << fabs(m.get_data(for_range+1,y_ind))/fabs(avg) << "\t" << plateau_condition << "\t" << for_range << endl;
        for_range++;
        avg = 0;//cin >> avg;
    }
    for_range--;
    ///for_range - ������ �����
    ///cout << "PLATEAU of " << max_plateau_condition << " level FROM " << for_range << "\t = " << m.get_data(for_range,y_ind) << endl;
    delete [] arr_of_avg;
    return for_range;
};

///FILTERS
matrix Get_Envelope_both1(matrix A, matrix dA)
{
    int arr[dA.get_height()]; ///������� ��������� �����, "�����"
    int j = 0;
    for(int i = 1; i < dA.get_height()-1; i++)
    {
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

    matrix E(j, 2);
    for(int i = 0; i < E.get_height()-1; i++)
    {
        E.set_data(i, 0, A.get_data(arr[i], 0));
        E.set_data(i, 1, A.get_data(arr[i], 1));
    }
    return E;
};

matrix Get_Envelope_top(matrix A, matrix dA)
{
    int arr[dA.get_height()];/// ������� ��������� �����, "�����"
    int j = 0;
    for(int i = 1; i < dA.get_height(); i++)
        if(dA.get_data(i-1,1) > 0.0) ///������ ���������
            if(dA.get_data(i,1) <= 0.0) ///������ ���������
            {
                arr[j] = i;//cout << j << "\t" << i << "\t" << dA.get_data(i,0) << "\t" << dA.get_data(i,1) << endl;
                j++;
            }//    cout << "j = " << j << endl;
    matrix E(j, 2);//E.info();//cout << arr[E.get_height()-1] << endl;
    for(int i = 0; i < E.get_height(); i++)
    {
        if(i!=0 && i!= E.get_height()-1)
        {
            if(abs(A.get_data(arr[i]-1, 1)) >= abs(A.get_data(arr[i], 1)))
                arr[i] -=1;
            if(abs(A.get_data(arr[i]+1, 1)) > abs(A.get_data(arr[i], 1)))
                arr[i] +=1;
                E.set_data(i, 0, A.get_data(arr[i], 0));
                E.set_data(i, 1, A.get_data(arr[i], 1));
        }else{
            E.set_data(i, 0, A.get_data(arr[i], 0));
            E.set_data(i, 1, A.get_data(arr[i], 1));
        }
    }
    return E;
};

matrix Get_Envelope_top_hard(matrix A, matrix dA, double level)
{
    int arr[dA.get_height()];/// ������� ��������� �����, "�����"
    int j = 0;
    for(int i = 1; i < dA.get_height(); i++)
        if(dA.get_data(i-1,1) > 0.0) ///������ ���������
            if(dA.get_data(i,1) <= 0.0) ///������ ���������
                if(A.get_data(i, 1) >= level)
                {
                    arr[j] = i;//cout << j << "\t" << i << "\t" << dA.get_data(i,0) << "\t" << dA.get_data(i,1) << endl;
                    j++;
                }//    cout << "j = " << j << endl;
    matrix E(j, 2);//E.info();//cout << arr[E.get_height()-1] << endl;
    for(int i = 0; i < E.get_height(); i++)
    {
        if(i!=0 && i!= E.get_height()-1)
        {
            if(abs(A.get_data(arr[i]-1, 1)) >= abs(A.get_data(arr[i], 1)))
                arr[i] -=1;
            if(abs(A.get_data(arr[i]+1, 1)) > abs(A.get_data(arr[i], 1)))
                arr[i] +=1;
                E.set_data(i, 0, A.get_data(arr[i], 0));
                E.set_data(i, 1, A.get_data(arr[i], 1));
        }else{
            E.set_data(i, 0, A.get_data(arr[i], 0));
            E.set_data(i, 1, A.get_data(arr[i], 1));
        }
    }
    return E;
};

matrix Get_Envelope_top(matrix A, int index_y, int index_x, matrix dA)
{
    int arr[dA.get_height()];/// ������� ��������� �����, "�����"
    int j = 0;
    for(int i = 1; i < dA.get_height(); i++)
        if(dA.get_data(i-1,1) > 0.0) ///������ ���������
            if(dA.get_data(i,1) <= 0.0)
            {
                arr[j] = i;//cout << j << "\t" << i << "\t" << dA.get_data(i,0) << "\t" << dA.get_data(i,1) << endl;
                j++;
            }//cout << "j = " << j << endl;

    matrix E(j, 2);//E.info();//cout << arr[E.get_height()-1] << endl;
    for(int i = 0; i < E.get_height(); i++)
    {
        if(i!=0 && i != E.get_height()-1)
        {
            if(abs(A.get_data(arr[i]-1, index_y)) >= abs(A.get_data(arr[i], index_y)))
                arr[i] -=1;
            if(abs(A.get_data(arr[i]+1, index_y)) > abs(A.get_data(arr[i], index_y)))
                arr[i] +=1;
                E.set_data(i, 0, A.get_data(arr[i], index_x));
                E.set_data(i, 1, A.get_data(arr[i], index_y));
        }else{
            E.set_data(i, 0, A.get_data(arr[i], index_x));
            E.set_data(i, 1, A.get_data(arr[i], index_y));
        }
    }
    return E;
}


matrix Get_Envelope_top_hard(matrix A, int index_y, int index_x, matrix dA, double level)///level - ������� ������ �������� ����������� ��������� ��������� (�������� ������� ������ �������� ����������� ���������)
{
    int arr[dA.get_height()];/// ������� ��������� �����, "�����"
    int j = 0;
    for(int i = 1; i < dA.get_height(); i++)
        if(dA.get_data(i-1,1) > 0.0) ///������ ���������
            if(dA.get_data(i,1) <= 0.0)
                if(A.get_data(i, index_y) >= level)
                {
                    arr[j] = i;//cout << j << "\t" << i << "\t" << dA.get_data(i,0) << "\t" << dA.get_data(i,1) << endl;
                    j++;
                }//    cout << "j = " << j << endl;

    matrix E(j, 2);//E.info();//cout << arr[E.get_height()-1] << endl;
    for(int i = 0; i < E.get_height(); i++)
    {
        if(i!=0 && i != E.get_height()-1)
        {
            if(abs(A.get_data(arr[i]-1, index_y)) >= abs(A.get_data(arr[i], index_y)))
                arr[i] -=1;
            if(abs(A.get_data(arr[i]+1, index_y)) > abs(A.get_data(arr[i], index_y)))
                arr[i] +=1;
                E.set_data(i, 0, A.get_data(arr[i], index_x));
                E.set_data(i, 1, A.get_data(arr[i], index_y));
        }else{
            E.set_data(i, 0, A.get_data(arr[i], index_x));
            E.set_data(i, 1, A.get_data(arr[i], index_y));
        }
    }
    return E;
}

matrix Get_Envelope_top(matrix y_column, matrix x_column, matrix dy_column)/// f(x)//x//df(x)/dx ///������������� ��� ������������ ������ ����!!!
{cout << "Get_Envelope_top(matrix y_column, matrix x_column, matrix dy_column)\n";
    //y_column.info(); x_column.info();dy_column.info();
    int arr[dy_column.get_height()];/// ������� ��������� �����
    int j = 0;//cout << dy_column << endl;
    for(int i = 1; i < dy_column.get_height(); i++)///����������� ����� �������� ����� ���� �����������
        if(dy_column.get_data(i-1,0) >= 0.0) ///���������� ���������
            if(dy_column.get_data(i,0) < 0.0 && y_column.get_data(i,0)>=0)
            {
                arr[j] = i;//cout << j << "\t"; cout << i << "\t" << y_column.get_data(i,0) << "\t" << dy_column.get_data(i,0) << endl;
                j++;
            }//    cout << "j = " << j << endl;

    ///��������� ������� ��� ����������� ��������� ����� 0, ������� �������.
    matrix E(j, 2);//E.info();//cout << arr[E.get_height()-1] << endl;
    for(int i = 0; i < E.get_height(); i++)
    {
        if(i!=0 && i != E.get_height()-1)
        {   //cout << arr[i] << "\t";
            arr[i]-=1;
            if(y_column.get_data(arr[i]+1,0) > y_column.get_data(arr[i],0))
            { //cout << x_column.get_data(arr[i],0) << "\t" << y_column.get_data(arr[i]+1,0) << "\t" << y_column.get_data(arr[i],0) << endl;
                arr[i]+=1;}
                E.set_data(i, 0, x_column.get_data(arr[i], 0));
                E.set_data(i, 1, y_column.get_data(arr[i], 0));
        }else{
            E.set_data(i, 0, x_column.get_data(arr[i], 0));
            E.set_data(i, 1, y_column.get_data(arr[i], 0));
        }
    }
    return E;
}

matrix Get_Envelop_frq_top(matrix x, matrix y, double freq, double ambit)
{cout << "matrix Get_Envelop_frq_top(matrix x, matrix y, double freq, double ambit)\n";
    double T = 1.0/freq; ///��� �� ������ � ���������� ��� ���������
       // T *= 1.0+ambit;

    matrix env;
    matrix buff(1,2);///����� ���������
    ///    env = env.merge_height(buff, env.get_height());

    ///step 1. Check before T
    int index_from = 0, index_to = 0, i = 0;
    while(x.get_data(index_to,0) < x.get_data(i,0) + T *(1.0+ambit))///� ���������� �� ������ �� �������
        index_to++;
    i = y.get_max_elind_column_gap(0, index_from, index_to); ///����� ������ ���������
        buff.set_data(0,0, x.get_data(i,0));
        buff.set_data(0,1, y.get_data(i,0));
            env = env.merge_height(buff, env.get_height());///��������� ��������� �����

while(x.get_data(i,0) + T *(1.0+ambit) < x.get_data(x.get_height()-1, 0))
{//cout << "2st stage\n";
    ///step Cycle. Check around T step
    index_from = i;
    while(x.get_data(index_from,0) < x.get_data(i,0) + T *(1.0-ambit))///����� ����� �����������
        index_from++;
           // cout << index_from << "\t";

    index_to = index_from;
    while(x.get_data(index_to,0) < x.get_data(i,0) + T *(1.0+ambit))///������ ����� �����������
        index_to++;
           // cout << index_to << "\t";

    i = y.get_max_elind_column_gap(0, index_from, index_to); //cout << i << "\t";
        buff.set_data(0,0, x.get_data(i,0));
        buff.set_data(0,1, y.get_data(i,0));
    env = env.merge_height(buff, env.get_height());///��������� ��������� �����
    //cout << buff << "/n";
}//cout << "end\n";
return env;
};

matrix Get_Envelop_frq_bot(matrix x, matrix y, double freq, double ambit)
{cout << "matrix Get_Envelop_frq_bot(matrix x, matrix y, double freq, double ambit)\n";
    double T = 1.0/freq; ///��� �� ������ � ���������� ��� ���������
       // T *= 1.0+ambit;

    matrix env;
    matrix buff(1,2);///����� ���������
    ///    env = env.merge_height(buff, env.get_height());

    ///step 1. Check before T
    int index_from = 0, index_to = 0, i = 0;
    while(x.get_data(index_to,0) < x.get_data(i,0) + T *(1.0+ambit))
        index_to++;
    i = y.get_min_elind_column_gap(0, 0, index_to); ///����� ��������
        buff.set_data(0,0, x.get_data(i,0));
        buff.set_data(0,1, y.get_data(i,0));
            env = env.merge_height(buff, env.get_height());///��������� ��������� �����
           // cout << env << "\t" << T << "\n1st stage\n";
    //cout << x.get_data(i,0) + T *(1.0+ambit) << "\t" << x.get_data(x.get_height()-1, 0) << endl;

while(x.get_data(i,0) + T *(1.0+ambit) < x.get_data(x.get_height()-1, 0))
{//cout << "2st stage\n";
    ///step 2. Check around T step
    index_from = i;
    while(x.get_data(index_from,0) < x.get_data(i,0) + T *(1.0-ambit))///����� ����� �����������
        index_from++;

    index_to = index_from;
    while(x.get_data(index_to,0) < x.get_data(i,0) + T *(1.0+ambit))///������ ����� �����������
        index_to++;

    i = y.get_min_elind_column_gap(0, index_from, index_to); //cout << i << "\t";
        buff.set_data(0,0, x.get_data(i,0));
        buff.set_data(0,1, y.get_data(i,0));
    env = env.merge_height(buff, env.get_height());///��������� ��������� �����
}//cout << "end\n";
return env;
};


matrix Get_Envelope_bot(matrix y_column, matrix x_column, matrix dy_column)/// f(x)//x//df(x)/dx ///������������� ��� ������������ ������ ����!!!
{//cout << "Get_Envelope_top(matrix y_column, matrix x_column, matrix dy_column)\n";
    y_column.info();
    x_column.info();
    dy_column.info();
    int arr[dy_column.get_height()];/// ������� ��������� �����
    int j = 0;
    for(int i = 1; i < dy_column.get_height(); i++)
        if(dy_column.get_data(i-1,0) < 0.0 && y_column.get_data(i,0)<=0) ///������ ���������
            if(dy_column.get_data(i,0) >= 0.0 )
            {
                arr[j] = i;//cout << j << "\t" << i << "\t" << y_column.get_data(i,0) << "\t" << dy_column.get_data(i,1) << endl;
                j++;
            }//    cout << "j = " << j << endl;

    matrix E(j, 2);//E.info();//cout << arr[E.get_height()-1] << endl;
    for(int i = 0; i < E.get_height(); i++)
    {
        if(i!=0 && i != E.get_height()-1)
        {
            arr[i]-=1;
            if(y_column.get_data(arr[i]+1,0) < y_column.get_data(arr[i],0))
            { //cout << x_column.get_data(arr[i],0) << "\t" << y_column.get_data(arr[i]+1,0) << "\t" << y_column.get_data(arr[i],0) << endl;
                arr[i]+=1;
            }

                E.set_data(i, 0, x_column.get_data(arr[i], 0));
                E.set_data(i, 1, y_column.get_data(arr[i], 0));
        }else{
            E.set_data(i, 0, x_column.get_data(arr[i], 0));
            E.set_data(i, 1, y_column.get_data(arr[i], 0));
        }
    }
    return E;
}

matrix Get_Envelope_top_hard(matrix y_column, matrix x_column, matrix dy_column, double level)/// f(x)//x//df(x)/dx
{
    int arr[dy_column.get_height()];/// ������� ��������� �����
    int j = 0;
    for(int i = 1; i < dy_column.get_height(); i++)
        if(dy_column.get_data(i-1,0) > 0.0) ///������ ���������
            if(dy_column.get_data(i,0) <= 0.0)
                if(y_column.get_data(i, 0) >= level)
                {
                    arr[j] = i;//cout << j << "\t" << i << "\t" << dA.get_data(i,0) << "\t" << dA.get_data(i,1) << endl;
                    j++;
                }//    cout << "j = " << j << endl;

    matrix E(j, 2);//E.info();//cout << arr[E.get_height()-1] << endl;
    for(int i = 0; i < E.get_height(); i++)
    {
        if(i!=0 && i != E.get_height()-1)
        {
            if(abs(y_column.get_data(arr[i]-1, 0)) >= abs(y_column.get_data(arr[i], 0)))
                arr[i] -=1;
            if(abs(y_column.get_data(arr[i]+1, 0)) > abs(y_column.get_data(arr[i], 0)))
                arr[i] +=1;
                E.set_data(i, 0, x_column.get_data(arr[i], 0));
                E.set_data(i, 1, y_column.get_data(arr[i], 0));
        }else{
            E.set_data(i, 0, x_column.get_data(arr[i], 0));
            E.set_data(i, 1, y_column.get_data(arr[i], 0));
        }
    }
    return E;
}


matrix Get_Envelope_bot(matrix A, matrix dA)
{
int arr[dA.get_height()];/// ������� ��������� �����

    int j = 0;
    for(int i = 1; i < dA.get_height()-1; i++)
        if(dA.get_data(i-1,1) <= 0.0)
            if(dA.get_data(i,1) >= 0.0)
            {
                arr[j] = i;
                j++;
            }

    matrix E(j, 2);
    for(int i = 0; i < E.get_height(); i++)
    {// cout << arr[i] << "\t" << this->angle.get_data(arr[i], 1) << "\t" << this->angle.get_data(arr[i], 0) << endl;
        if(i!=0 && i!= E.get_height()-1)
        {
            if(abs(A.get_data(arr[i]-1, 1)) >= abs(A.get_data(arr[i], 1)))
                arr[i] -=1;
            if(abs(A.get_data(arr[i]+1, 1)) > abs(A.get_data(arr[i], 1)))
                arr[i] +=1;
                E.set_data(i, 0, A.get_data(arr[i], 0));
                E.set_data(i, 1, A.get_data(arr[i], 1));
        }else{
            E.set_data(i, 0, A.get_data(arr[i], 0));
            E.set_data(i, 1, A.get_data(arr[i], 1));
        }
    }
    return E;
};


matrix Get_Envelope_bot(matrix A, int index_y, int index_x, matrix dA)
{
int arr[dA.get_height()];/// ������� ��������� �����

    int j = 0;
    for(int i = 1; i < dA.get_height()-1; i++)
        if(dA.get_data(i-1,1) <= 0.0)
            if(dA.get_data(i,1) >= 0.0)
            {
                arr[j] = i;
                j++;
            }

    matrix E(j, 2);
    for(int i = 0; i < E.get_height(); i++)
    {// cout << arr[i] << "\t" << this->angle.get_data(arr[i], 1) << "\t" << this->angle.get_data(arr[i], 0) << endl;
        if(i!=0 && i!= E.get_height()-1)
        {
            if(abs(A.get_data(arr[i]-1, index_y)) >= abs(A.get_data(arr[i], index_y)))
                arr[i] -=1;
            if(abs(A.get_data(arr[i]+1, index_y)) > abs(A.get_data(arr[i], index_y)))
                arr[i] +=1;
                E.set_data(i, 0, A.get_data(arr[i], index_x));
                E.set_data(i, 1, A.get_data(arr[i], index_y));
        }else{
            E.set_data(i, 0, A.get_data(arr[i], index_x));
            E.set_data(i, 1, A.get_data(arr[i], index_y));
        }
    }
    return E;
};

matrix Get_Envelope_topbot(int index_govern, matrix env_t, matrix env_b)
{cout << "Get_Envelope_topbot(int index_govern, matrix env_t, matrix env_b)\n";
    matrix result(env_t.get_height() + env_b.get_height()-1, env_t.get_width());
    for(int i = 0, j = 0, k = 0; /*i < env_t.get_height(), j< env_b.get_height(),*/ k < result.get_height();)
        if(env_t.get_data(i, index_govern)< env_b.get_data(j, index_govern))
        {
            for(int w = 0; w < env_t.get_width(); w++)
                result.set_data(k, w, env_t.get_data(i,w));
            k++;i++;
        }else{
            for(int w = 0; w < env_b.get_width(); w++)
                if(w == index_govern)
                    result.set_data(k, w, env_b.get_data(j,w));
                else
                    result.set_data(k, w, env_b.get_data(j,w)*(-1));///������ ������������� ������������ ����
            k++;j++;
        }
    return result;
}

void Window_Sigma_Filter(double* arr, const int arr_length, int window_range, const double sigma_level)
{cout << "Window_Sigma_Filter(double* arr, const int arr_length, int window_range, const int ejection_level)\n";
    cout << "window_range (>)= " << window_range << "\tejection_level = " << sigma_level << endl;
    if(window_range < arr_length)
    {
    double sum = 0, avg = 0, sigma = 0;
        sum = window_range;

    if(window_range/2 == sum/2)///���� ������, �� �������� �� ��������� � �������� � ���
            window_range++;
            //cout << "window_range = " << window_range << endl;
    sum = 0;

    double n_arr [arr_length];
    bool fixing_map [arr_length];
    ///window_range ��������� ��� ������/��������
        for(int j = (window_range-1)/2; j < arr_length - (window_range-1)/2; j++)///�������� ��� ����� ����� ���� ������� ������� ///+1 �� 2 ��������� ��� ����������� �������� ����� ���� ���������
        {//cout << "j = " << j << endl;
            ///�������� � ����, �������� ��������������� �����
            for(int i = j - (window_range-1)/2; i < j + ((window_range-1)/2); i++)
                if(i!=j){sum += arr[i];}
            avg = sum/(window_range-1);

            ///����� � ����� j = �������� ����
            sigma = sqrt(pow(arr[j] - avg, 2)/(window_range-1));
            //cout << arr[j - (window_range-1)/2]  << "\t" << arr[j] << "\t" << arr[j +((window_range-1)/2)] << "\t sigma = " << sigma << endl;

            if (sigma != 0) cout << j << "\t" << sigma << endl;
            ///������ �����
            if( fabs(fabs(arr[j]) - fabs(avg)) >= sigma_level*sigma && sigma != 0)/// true ���� ����� ������ //if( sigma >= ejection_level)
            {
                fixing_map[j] = true;
                    cout << "fix in " << j << "\tsigma\t" << sigma << endl;
            }else{
                fixing_map[j] = false;
            }
               // cout << "fixing_map = " << fixing_map[j] << endl;
            ///������� � ���������� ���� �����������
            sum = 0;
        }
    for(int i = 0 ; i < arr_length; i++)
        n_arr[i] = arr[i];
    for(int i = 0 ; i < arr_length; i++)
    {
        if(fixing_map[i])
        {
            arr[i] = (n_arr[i-1] + n_arr[i+1])/2.0; cout << "fixing " << i << "\t" << arr[i] << endl;
        }
    }
    }else{
    cout << "window_range >= arr_length ! procedure aborted\n";
    }
};



void Filter_Outlier1(double* in_arr, int from, int arr_length, double in_coeff, double out_coeff)
{//cout << "arr_length = " << arr_length << endl;

    double arr [arr_length];
    for(int i = 0; i < arr_length; i++)
        arr[i] = in_arr[i];

    bubble_sort(arr, arr_length);

    int num1 = 0, num2 = 0, num3 = 0;

    ///Q2
    double Q2 = 0; ///������� �������
    if(arr_length%2 == 0)
    {
        ///������
        num2 = arr_length/2-1; ///�� ����� ��������� � �����
        Q2 = (arr[num2] + arr[num2+1])/2; /// ����� ��������� � ��������
    }else{
        num2 = (int)arr_length/2; ///�������� � ����� �������� ����� � �������� // �� ��������� � ����
        Q2 = arr[num2]; ///��������
    }//cout << "Q2 = " << Q2 << "\t" << num2 << endl;

    double Q1 = 0; ///������� �������
    if((arr_length - num2)%2 == 0)///� ������ ��������
    {
       num1 = (int)(arr_length - num2)/2-1;
       Q1 = (arr[num1] + arr[num1+1])/2; ///������
    }else{
        num1 = (int)(arr_length - num2)/2;
        Q1 = arr[num1]; ///��������
    }//cout << "Q1 = " << Q1 << "\t" << num1 << endl;

    ///quartile
    double Q3 = 0;
    if((arr_length-1 - num2)%2 == 0)
    {//cout << "even\t" << (double)(arr_length-1 - num2)/2 << endl;
        ///������
        num3 = (int)arr_length/2-1 + int(num2/2+1); ///�� ������ ��������
        Q3 = (arr[num3] + arr[num3+1])/2;
    }else{//cout << (double)(arr_length-1 - num2)/2 << "\t" << (int) (arr_length-1 - num2)/2 << endl;
        ///��������
        num3 = (int)arr_length/2 + int(num2/2+1);//num2;
        Q3 = arr[num3];
    }//cout << "Q3 = " << Q3 << "\t" << num3 << endl;

    double QQ = Q3 - Q1;

    double in_gap1 = 0, in_gap2 = 0;
    in_gap1 = Q1  - (QQ)*in_coeff;
    in_gap2 = Q3  + (QQ)*in_coeff;

    double out_gap1 = 0, out_gap2 = 0;//cout << "out coeff = " << out_coeff << endl;
    out_gap1 = Q1 - (QQ)*out_coeff;
    out_gap2 = Q3 + (QQ)*out_coeff;

    double buff_arr[arr_length];
        buff_arr[0] = in_arr[0];
        buff_arr[arr_length-1] = in_arr[arr_length-1];

    ///������ ����� ///���������� ������
    if(in_arr[0] <= out_gap1 || in_arr[0] >= out_gap2)
        buff_arr[0] = in_arr[1];
    for(int i = 1; i < arr_length-1; i++) ///��� ���� �������
    {//cout << in_arr[i] << "\t" << out_gap1 << "\t" << out_gap2 << endl;
         buff_arr[i] = in_arr[i];
        if(in_arr[i] <= out_gap1 || in_arr[i] >= out_gap2) ///���� ����� �������� ������ �� ���� ������ �� ������ �������
        {
            buff_arr[i] = (in_arr[i-1] + in_arr[i+1])/2;
            ///in_arr[i] = (in_arr[i-1] + in_arr[i+1])/2;
            //cout << i << "\tFILTERED\t" << buff_arr[i] << endl;
        }
        //arr[i-1] = buff1;
    }
    ///��������� ����� ///���������� �������������
    if(in_arr[arr_length-1] <= out_gap1 || in_arr[arr_length-1] >= out_gap2)
        {buff_arr[arr_length-1] = in_arr[arr_length-2];cout << arr_length-1 << "\t" <<  arr_length << endl;}

    for(int i = 0; i < arr_length; i++)
        in_arr[i] = buff_arr[i];
};

void Filter_Outlier1(double* in_arr, int arr_length, double in_coeff, double out_coeff)
{
    Filter_Outlier1(in_arr, 0, arr_length, in_coeff, out_coeff);
}

void Filter_Outlier(double* in_arr, int arr_length, int window_range, double in_coeff, double out_coeff) ///���������� ��� �� ������
{
    int i = 0;
    for( ; i < arr_length;)
    {//cout << "i =" << i << endl;
        if(arr_length - i < window_range)
        {   //cout << "IFFFFFFFFFFFFFFFF\n";
            ///����� �����
            in_arr -= (window_range - (arr_length - i));
            i = arr_length - window_range;
            break;
        }

        Filter_Outlier1(in_arr, 0, window_range, in_coeff, out_coeff);

        in_arr += window_range;
        i +=window_range;
    }
}

void smooth(double *input, double *output, int n, int window)
{
   int i,j,z,k1,k2,hw;
   double tmp;
   if(fmod(window,2)==0) window++;
   hw=(window-1)/2;
   output[0]=input[0];

   for (i=1;i<n;i++){
       tmp=0;
       if(i<hw){
           k1=0;
           k2=2*i;
           z=k2+1;
       }
       else if((i+hw)>(n-1)){
           k1=i-n+i+1;
           k2=n-1;
           z=k2-k1+1;
       }
       else{
           k1=i-hw;
           k2=i+hw;
           z=window;
       }

       for (j=k1;j<=k2;j++){
           tmp=tmp+input[j];
       }
       output[i]=tmp/z;
   }
}

matrix smooth_matrix_column(matrix A, int window_length, int column_index)
{
    if(window_length>1)
    {
        double arr[A.get_height()], n_arr[A.get_height()];
        ///������
        for(int i = 0; i < A.get_height(); i++)
            arr[i] = A.get_data(i, column_index);
    matrix n_A = A;
    smooth(arr, n_arr, A.get_height(), window_length);
        n_A.insert_array_to_matrix_column(column_index, n_arr, A.get_height());//cout << A.merge_width(n_A, A.get_width()) << endl;
    return n_A;
    }
    cout<< "smooth_matrix_column(matrix A, int window_length, int column_index) WARNING same data returned\n";
    return A;
}

matrix smooth_matrix_column1(matrix A, int window_length, int column_index)
{
    if(window_length>1)
    {
        double arr[A.get_height()], n_arr[A.get_height()];
        ///������
        for(int i = 0; i < A.get_height(); i++)
            arr[i] = A.get_data(i, column_index);
    matrix n_A = A;
    smooth(arr, n_arr, A.get_height(), window_length);
        n_A.insert_array_to_matrix_column(column_index, n_arr, A.get_height());//cout << A.merge_width(n_A, A.get_width()) << endl;
    return n_A;
    }
    cout<< "smooth_matrix_column(matrix A, int window_length, int column_index) WARNING same data returned\n";
    return A;
}

matrix smooth_matrix_height(matrix A, int window_length)
{
    for(int i = 0; i < A.get_width(); i++)
        A = smooth_matrix_column(A, window_length, i);
    return A;
}
///****

matrix smooth_matrix_height(matrix A, matrix window_lengthes)
{
    for(int i = 0; i < A.get_width(); i++)
        A = smooth_matrix_column(A, window_lengthes.get_data(0, i), i);//cout << "smooth_matrix_height\t" << i << endl;
    return A;
}

double f_poly(double f_arg, double* coeff, int power)///����������� � ������� ���������
{
    double result = 0;
    while(power >= 0)
    {
        result += pow(f_arg, power)*coeff[power];
            power--;
    }
    return result;
}

double Div_Dif(double* f_arr, int f_arr_length, double* arg_arr, int arg_arr_length, double* dd_arr, int& dd_arr_length)
{
    int N = arg_arr_length + 1;
    ///1 + 2 + 4 + ... + 2*N-1 = 1 + summ(2, N) an=a1+(n-1)d; a1 = 0; d=2
    ///����� ��������� � ������������
    dd_arr_length = 1 + (2*2+2*(N-1))/2*N;
    cout << dd_arr_length << endl;
    dd_arr = new double [dd_arr_length];
    int dd_counter = 0;
    for(int i = 0; i < 2*N - 1; i++)
    {
        for(int j = 0; j < N; j++)
        {
            dd_arr[dd_counter] = f_poly(arg_arr[i], f_arr, f_arr_length);
            dd_counter++;
        }
    }
}

void Newton_Polinom(double* arr, int arr_length, double* coeff, int polynom_power)
{
    coeff = new double [polynom_power];///0, 1 ... polynom_power
    double* dd = NULL;
    int dd_length = 0;

    Div_Dif(coeff, polynom_power, arr, arr_length, dd, dd_length);
};

double Polynom_Lagrange(double argument, double* arr_points, int poly_index, int points_amount)
{//cout << "Polynom_Lagrange\n";
    //cout << "poly_index " << poly_index << "\tpoly_arg " << argument << "\t";
    double r = 1.0;
    for(int i = 0; i < poly_index; i++)
    {
        r *= (argument - arr_points[i]);//cout << "l = " <<  r << "\t" << "l = " << i << "\t";
    }
    for(int i = poly_index + 1; i < points_amount; i++)
    {
        r *= (argument - arr_points[i]);//cout << "l = " <<  r << "\t" << "l = " << i << "\t";
    }//cout << "r = " << r << "\t";
    return r;
};


double Get_Lagrange(double argument, double* arr_func, double* arr_points, int points_amount)
{
    double r = 0;
    double P_a, P_i;
    //for(int i = 0; i < points_amount; i++)
       // cout << i << "\t" << arr_func[i] << endl << endl;
    for(int i = 0; i < points_amount; i++)
    {//cout << i << endl;
        P_a = Polynom_Lagrange(argument, arr_points, i, points_amount); //cout << P_a << "\t";
        P_i = Polynom_Lagrange(arr_points[i], arr_points, i, points_amount); //cout << P_i << "\n";
        if(P_i != 0)
            r += (arr_func[i] * P_a/P_i);//cout << endl;
    }
    return r;
};


/*double Approximatin(double argument, double* y_data, double* x_data, int approx_data_length, double (*approx_func)(double))
{
    ///case
    ///1, 2, 3 - poly power
    ///sin, cos, tg, ctg ....
    ///#include "matrix.h"
    ///approximate from data from 0 to approx_data_length
    return 0;
};

void Approximation_Window(double* y,  double* x, double* argument, int data_length, int window_length, double (*approx_func)(double))///������� ������������� �� ����� ����������
{///y - ��������� �������� ������� ��� ��������� x, argument - ����� � ������ ���������� ����� �������� �������

    if(window_length <= data_length)///��������������� �������
    {
        ///if(window_length == data_length)///������������� �� ����� ����������
        ///������ ��� �������������
        int n_data_length = window_length/2 * 2 + 1;
        double n_y [n_data_length];
        double n_x [n_data_length];
        ///���������� ������ ��� �������������
        for(int i = 0; i < window_length/2; i++)
        {
            n_y[i] = y[0];
            n_x[i] = x[0];
        }
        for(int i = window_length/2; i < n_data_length - window_length/2; i++)
        {
            n_y[i] = y[i];
            n_x[i] = x[i];
        }
        for(int i = n_data_length - window_length/2; i < n_data_length; i++)
        {
            n_y[i] = y[data_length-1];
            n_x[i] = x[data_length-1];
        }
        ///��������� ������� �� �������������� ������ �� ���� window_length
        for(int i = 0; i < data_length; i++)
        {
            ///��������� ������ ������� ������ ����������
            ///poly?
            ///������� ������� ������� - ������� ������
            y[i] = Approximatin(argument[i], n_y + i, n_x + i, window_length, approx_func);///����� ������� ���������� ����������� ������� ��� �������������, �����
        }
    }
}*/


double Standart_Deviation(double* arr, double* sigma, int data_length,  int window_length)
{
    double avg = 0;
    double sum1 = 0;
    for(int i = 0; i < data_length; i++)
    {
        for(int j = i; j < window_length + i; j++ )
            sum1 = pow(arr[i] - avg, 2);
        sigma[i] = sqrt(sum1/(data_length-1));
    }
};

/*
double Student_Coeff(double P, int f)
{
    if(P = 0.90)
        return 2;
    cout << "WARNING double Student_Coeff(double P, int f) 1 returned\n";
    return 1;
}

double part_dispertion(double* arr, int data_length)
{
    double avg = 0.0;
    for(int i = 0; i < data_length; i++)
        avg += arr[i];
    avg *= (1.0/data_length);

    double s = 0;
    double sum = 0;
    for(int i = 0; i < data_length; i++)
        sum = pow(arr[i]- avg, 2);
    s = sqrt(sum/(data_length - 1));
    return s;
    ///����������� + double P ������������� ��������, � ������� 0.90
    //int f = data_length - 1; ///Student coeff
    //double err = Student_Coeff(P, f)*s/sqrt(data_length);
};

double part_ref_dispertion(double* arr, int data_length)
{
    double avg = 0.0;
    for(int i = 0; i < data_length; i++)
        avg += arr[i];
    avg *= (1.0/data_length);

    return  part_dispertion(arr, data_length)/fabs(avg);
};

void arr_dispertion(double* arr, double* dispertion, int arr_length, int part_length)///���� �����������, ��������������
{
    double d;
    for(int i = 0; i < arr_length; )
    {
        if( i + part_length < arr_length)
        {
            d = part_dispertion(arr, part_length);
            for(int j = i; j < i + part_length; j++)
                dispertion[j] = d;
            ///������� � ��������� �����
            i += part_length;
            arr += part_length;
        }else{
            arr -= (i + part_length - arr_length);
            d = part_dispertion(arr, part_length);
            for(int j = i; j < arr_length; j++)
                dispertion[j] = d;
            i = arr_length;
        }
    }
}

void arr_dispertion(double* arr, double* dispertion, int arr_length, int window_range)///���� �����������, ��������������
{
    window_range = (window_range/2)*2+1;

    int d = 1;
    int i = 0;


    dispertion[i] = 0;
    d += 2;
    i++;

    for( ; i < window_range - (d-1)/2; i++) /// ������ �����, ����� ������ ������ � ����� j
    {
        dispertion[i] = part_dispertion(arr, d);
        d += 2;
    }

    for( ; i < arr_length - (window_range - (d-1)/2); i++)
    {
        dispertion[i] = part_dispertion(arr, d);
        arr++;
    }

    for( ; i < arr_length; i++)
    {
        dispertion[i] = part_dispertion(arr, d);
        d -= 2;
    }
}

void arr_ref_dispertion(double* arr, double* dispertion, int arr_length, int window_range)///���� �����������, ��������������
{
    window_range = (window_range/2)*2+1;

    int d = 1;
    int i = 0;


    dispertion[i] = 0;
    d += 2;
    i++;

    for( ; i < window_range - (d-1)/2; i++) /// ������ �����, ����� ������ ������ � ����� j
    {
        dispertion[i] = part_ref_dispertion(arr, d);
        d += 2;
    }

    for( ; i < arr_length - (window_range - (d-1)/2); i++)
    {
        dispertion[i] = part_ref_dispertion(arr, d);
        arr++;
    }

    for( ; i < arr_length; i++)
    {
        dispertion[i] = part_ref_dispertion(arr, d);
        d -= 2;
    }
}*/

matrix derevative_3_dot(matrix Y, int y_index, int x_index)
{///cout << "derevative_3_dot(matrix Y, int y_index, int x_index)" << endl;
    matrix D = Y.get_column(x_index);
    D = D.merge_width(Y.get_column(y_index), D.get_width());
    ///D.info();
    ///��� ����������� �����
    double h = Y.get_data(1,x_index) - Y.get_data(0,x_index);
    ///������
        D.set_data(0, 1, (-3*Y.get_data(0, y_index) + 4*Y.get_data(1, y_index) - Y.get_data(2, y_index))/2.0/h );
    ///�������������
    for(int i = 1; i < Y.get_height()-1; i++)
        D.set_data(i, 1, (Y.get_data(i+1, y_index) - Y.get_data(i-1, y_index))/2.0/h );//cout << "DLAST\n";
    ///���������
        D.set_data(D.get_height()-1, 1, (Y.get_data(Y.get_height()-3, y_index) - 4*Y.get_data(Y.get_height()-2, y_index) + 3*Y.get_data(Y.get_height()-1, y_index))/2.0/h);
    return D;
}

matrix derevative_3_dot_h(matrix Y, int y_index, int x_index)///�� aco.ifmo.ru/el_books/numerical_methods/lectures/glava1.html
{///cout << "derevative_3_dot(matrix Y, int y_index, int x_index)" << endl;
    matrix D = Y.get_column(x_index);
    D = D.merge_width(Y.get_column(y_index), D.get_width());
    ///D.info();
    double h = 0;
    ///������
        h = Y.get_data(2,x_index) - Y.get_data(0,x_index); ///����� ����� �� h ������  ///
        D.set_data(0, 1, (-3*Y.get_data(0, y_index) + 4*Y.get_data(1, y_index) - Y.get_data(2, y_index))/h );
    ///�������������
    for(int i = 1; i < Y.get_height()-1; i++)
    {
        h = Y.get_data(i+1,x_index) - Y.get_data(i-1,x_index);
        D.set_data(i, 1, (Y.get_data(i+1, y_index) - Y.get_data(i-1, y_index))/2.0/h );
    }//cout << "DLAST\n";
    ///���������
        h = Y.get_data(Y.get_height()-1,x_index) - Y.get_data(Y.get_height()-3,x_index);
        D.set_data(D.get_height()-1, 1, (Y.get_data(Y.get_height()-3, y_index) - 4*Y.get_data(Y.get_height()-2, y_index) + 3*Y.get_data(Y.get_height()-1, y_index))/2.0/h);
    return D;
}


matrix derevative_3_dot(matrix y, matrix x)
{//cout << "derevative_3_dot(matrix Y, int y_index, int x_index)" << endl;
    matrix D = x;
        D = D.merge_width(y, D.get_width());
    ///D.info();
    ///��� ����������� �����
    double h = D.get_data(1,0) - D.get_data(0,0);//cout << h << endl;
    ///������
        D.set_data(0, 1, (-3*y.get_data(0, 0) + 4*y.get_data(1, 0) - y.get_data(2, 0))/2.0/h );
    ///�������������
    for(int i = 1; i < D.get_height()-1; i++)
        D.set_data(i, 1, (y.get_data(i+1, 0) - y.get_data(i-1, 0))/2.0/h );
    ///���������
        D.set_data(D.get_height()-1, 1, (y.get_data(y.get_height()-3, 0) - 4*y.get_data(y.get_height()-2, 0) + 3*y.get_data(y.get_height()-1, 0))/2.0/h);
    return D;
}

double factorial2(double arg)
{
double r = 1;
    if(arg == 0)
        return 1;
    if(((int)arg/2)*2 == arg)
    {
        while(arg>0)
        {
            r *= arg;
            arg-=2;
        }
        return r;
    }else{
        while(arg>1)
        {
            r *= arg;
            arg-=2;
        }
        return r;
    }
return r;
};

double Legander_K(double arg, int degree)
{
    double K = 1.0, d = 0;
    for(int i = 2; i <= degree; )
    {
        d = pow(factorial2(i-1)/factorial2(i), 2)*pow(sin(arg), i);//cout <<  i << "\t" << d << "\t" << i-1 << "\t" << i << "\t" << factorial2(i-1)/factorial2(i) << endl;
            K = K + d;
        i+=2;
    }
    return K;
};

void equation_2order(const double a, const double b, const double c,
                        double& root1Re, double& root1Im, double& root2Re, double& root2Im)
{
    double D = b*b - 4.0*a*c;
    ///root1 +sqrt(D)
    if(D >= 0.0)
    {///onlyRe
        root1Re = (-1.0*b+sqrt(D))/2.0;
            root1Im = 0.0;
        root2Re = (-1.0*b-sqrt(D))/2.0;
            root2Im = 0.0;
    }else{///Re+Im
        root1Re = -1.0*b/2.0;
            root1Im = sqrt(-1.0*(D))/2.0;
        root2Re = -1.0*b/2.0;
            root2Im = -1.0*sqrt(-1.0*(D))/2.0;
    }
};


double matrix_func(matrix m, double arg, int index_arg, int index_func)///������� ������������ �������� - ��������� �������� ������� �� ������ �������
{
    for(int i = 0; i < m.get_height()-1; i++)
    {
        if(arg < m.get_data(i,index_arg))
            return m.get_data(i,index_func);
        if(arg > m.get_data(i,index_arg))
            if(arg < m.get_data(i+1,index_arg))
                return m.get_data(i,index_func);
    }
    if(arg > m.get_data( m.get_height()-1,index_arg))
        return m.get_data(m.get_height()-1,index_func);
}

int sign(double a)
{
    if(a < 0)return -1.0;
    if(a > 0) return 1.0;
    if(a == 0) return 0.0;
    return 0.0;
};
