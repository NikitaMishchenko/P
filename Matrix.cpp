#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <io.h>

#include "Matrix.h"

matrix::matrix() ///empty
{
    height = 0; width = 0;
    data = NULL;
};

void matrix::matrix_alloc(int h, int w)
{
    data = new double* [h];
    for (int i = 0; i < h; i++)
    {
        data[i] = new double[w];
            for (int j = 0; j < w; j++)
                data[i][j] = 0;
    }
};

matrix::matrix(int h, int w) ///rectangle
{
    height = h; width = w;
    this->matrix_alloc(h, w);
};

matrix::matrix(int n) ///square
{
    height = n; width = n;
    this->matrix_alloc(n, n);
};

matrix::matrix(int n, char c)
{
    height = n;
    width = n;
        this->matrix_alloc(n, n);
    if (c == 'E')
    {
        for(int i = 0; i < n; i++)
            data[i][i] = 1.0;
    }
    if (c == 'N')
    {
         for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                data[i][j] = i;
    }
};

matrix::matrix(int h, int w, char c)
{
    height = h;
    width = w;
        this->matrix_alloc(h, w);
    if (c == 'N')
    {
         for(int i = 0; i < h; i++)
            for(int j = 0; j < w; j++)
                data[i][j] = i;
    }
};

matrix::~matrix() ///destructor
{
    for (int i = 0; i < height; i++)
        delete [] data[i];
    delete [] data;
};

matrix::matrix(const matrix &copy_A) ///copy
{//cout << "copy_matrix\n";
    data = new double* [copy_A.height];
    for (int i = 0; i < copy_A.height; i++)
    {
        data[i] = new double[copy_A.width];
            for (int j = 0; j < copy_A.width; j++)
                data[i][j] = copy_A.data[i][j];
    }
    height = copy_A.height;
    width = copy_A.width;
}

///primary
void matrix::delete_data()
{
    for (int i = 0; i < this->height; i++)
        delete [] data[i];
            delete [] data;
    height = 0;width = 0;
    data = new double* [height];
    for (int i = 0; i < height; i++)
    {
        data[i] = new double[width];
    }
};

int matrix::get_height() const
    {return height;};

int matrix::get_width() const
    {return width;};

double matrix::get_data(int i, int j) const
    {return data[i][j];};

void matrix::set_height(int h)
    {height = h;};

void matrix::set_width(int w)
    {width = w;};

void matrix::set_data(int i, int j, double d)
    {data[i][j] = d;};

//void matrix::set_column_data(int h, int w, double data)
//{
//};

void matrix::make_matrix_data(int h, int w, double d)
{
    matrix res(h,w);
    res.info();
    for(int i =0; i < h; i++)
        for(int j =0; j < w; j++)
            res.set_data(i,j,d);
    *this = res;
};


void matrix::insert_array_to_matrix_column(int index_column, double* arr, int arr_length)///overwrite and fit to matrix.height();
{
    if(arr_length > this->get_height())
        std::cout << "matrix::insert_array_to_matrix_column(int index_column, double* arr, int length) WARNING array cut to matrix.height()\n";
    for(int i = 0; i < this->get_height(); i++)
        this->set_data(i,index_column, arr[i]);
}

void matrix::insert_array_to_matrix_row(int index_row, double* arr, int arr_length)
{
    if(arr_length > this->get_width())
        std::cout << "matrix::insert_array_to_matrix_column(int index_column, double* arr, int length) WARNING array cut to matrix.width()\n";
    for(int i = 0; i < this->get_width(); i++)
        this->set_data(index_row, i, arr[i]);
}

double matrix::get_max_el_column(int index_column)
{//cout << "double matrix::get_max_el_column(int index_column)\n";
if(height != 0)
    {double m = data[0][index_column];
    for(int i = 1; i < height; i++)
        if(data[i][index_column] > m)
            m = data[i][index_column];
    return m;}
return -999999;
};

double matrix::get_max_el_column_gap(int index_column, int gap_index_from, int gap_index_to)
{
if(height != 0)
    {double m = data[gap_index_from][index_column];
    for(int i = gap_index_from+1; i < gap_index_to; i++)
        if(data[i][index_column] > m)
            m = data[i][index_column];
    return m;}
return -999999;
}

double matrix::get_min_el_column(int index_column)
{//cout << "double matrix::get_min_el_column(int index_row)\n";
if(width != 0)
    {double m = data[0][index_column];
    for(int i = 1; i < height; i++)
        {
        if(data[i][index_column] < m)
            m = data[i][index_column];//cout << data[i][index_column] << endl;
        }
    return m;}
return -999999;
};


int matrix::get_max_elind_column_gap(int index_column, int gap_index_from, int gap_index_to)
{
    int index = gap_index_from;
if(height != 0)
    {double m = data[gap_index_from][index_column];
    for(int i = gap_index_from+1; i < gap_index_to-1; i++)
        if(data[i][index_column] > m)
        {
            m = data[i][index_column];
                index = i;
        }
    return index;}
return -999999;
};

int matrix::get_min_elind_column_gap(int index_column, int gap_index_from, int gap_index_to)
{
    int index = gap_index_from;
if(height != 0)
    {double m = data[gap_index_from][index_column];
    for(int i = gap_index_from+1; i < gap_index_to-1; i++)
        if(data[i][index_column] < m)
        {
            m = data[i][index_column];
                index = i;
        }
    return index;}
return -999999;
};

///OPERATORS
matrix& matrix::operator= (const matrix& A)
{//cout << "matrix::operator=\n";
if( this != &A)
{
    for (int i = 0; i < this->height; i++) ///changed width ->height
        delete [] data[i];
    delete [] data;

    height = A.height;
    width = A.width;
    data = new double* [height];
    for (int i = 0; i < height; i++)
    {
        data[i] = new double[width];
            for (int j = 0; j < width; j++)
                data[i][j] = A.data[i][j];
    }
}
return *this;
};

 matrix& matrix::operator+=(const matrix& A)
 {
    if(height== A.height && width == A.width)
    {
        for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    data[i][j] = data[i][j]+ A.data[i][j];
        return *this;
    }else{
        std::cout << "Warning! Different sizes of matrices << operator+= >> Same matrix returned" << std::endl;
        return *this;
    }
 };

 matrix operator+(const matrix& A, const matrix& B)
 {
    matrix result = A;
    result += B;
    return result;
 };

matrix matrix::operator* (double a)
{
    for(int i = 0; i < height; i++)
        for (int j = 0; j < width; j++)
            data[i][j] = a*data[i][j];
    return *this;
};

matrix operator*(const double& a, const matrix& A)
{
    matrix r(A.height, A.width);
    for(int i = 0; i < r.height; i++)
        for (int j = 0; j < r.width; j++)
            r.data[i][j] = a*A.data[i][j];
    return r;
}

matrix operator*(const matrix& A, const matrix& B)
{//cout << "operator * curruped\n";
///falls with no reason????? counted right
    matrix result(A.height, B.width);
if(A.width == B.height)
{
    double buff = 0.0;
    for(int i = 0; i < B.height; i++)
        for(int j = 0; j < B.width; j++)
        {
            buff = 0.0;
            for(int k = 0; k < B.height; k++)
                buff += A.data[i][k] * B.data[k][j];
            result.data[i][j] = buff;
        }
}
else{std::cout << "Wrong matrix sizes in operator* A.height = "
        << A.width  << "\tB.width = " << B.height << "\tEmpty matrix returned " << std::endl;}
    return result;
};

matrix matrix::operator*=(const matrix& A)
{
    if(height== A.width && width == A.height)
    {
        *this =  (*this * A);
        return *this;
    }else{
        std::cout << "Different sizes of matix in << operator*= >>" << std::endl;
        return *this;
    }
}

matrix matrix::operator/ (double a)
{
    return *this*(1.0/a);
};


matrix matrix::operator-= (const matrix& A)
{
    for (int i = 0; i < height; i++)
                for (int j = 0; j < width; j++)
                    data[i][j] = data[i][j]- A.data[i][j];
    return *this;
};

matrix operator-(const matrix& A, const matrix& B)
{
    matrix result = A;
        result -= B;
    return result;
};

std::ostream & operator << (std::ostream &C, const matrix &A)
{
    for (int i = 0; i < A.height; i++)
    {
        for (int j = 0; j < A.width; j++)
            C << A.data[i][j] << "\t";
        if( i!= A.height-1 )
            C << "\n";
    }
    return C;
 };

std::istream & operator >> (std::istream &C, matrix &A)
{
    for (int i = 0; i < A.height; i++)
        for (int j = 0; j < A.width; j++)
            C >> A.data[i][j];
    return(C);
};


matrix matrix::abs_column(int index_column)
{
    matrix res = *this;
    for(int i = 0; i < height; i++)
        res.set_data(i,index_column, fabs(this->get_data(i,index_column)));
    return res;
};

void matrix::abs_index_column(int index_column)
{
    for(int i = 0; i < height; i++)
        this->set_data(i,index_column, fabs(this->get_data(i,index_column)));
}

matrix matrix::transpose()
{
    matrix T(this->width, this->height);
        for (int i = 0;  i < this->height; i++)
            for(int j = 0; j < this->width; j++)
                T.data[j][i] = this->data[i][j];
    return(T);
};

void matrix::transpose(matrix& m)
{
    *this = m.transpose();
};

matrix matrix::minor(int l,int k)
{
    int n = this->height-1;
    int s1 = 0; int s2 = 0; int p = 0;
    matrix M(n);
    for (int i = 0; i < n; i++)
    {
        if(i==l){s1=1;};
        p=0;
        for (int j = 0; j < n; j++)
        {
            if (j==k){s2=1; p++;};
            if(p%2 == 0 ){s2=0;};
            M.data[i][j]=this->data[i+s1][j+s2];
        }
    }
    return(M);
};

double matrix::determinant()
{
    double d = 0;
if(this->height == this->width)
{
    int n = this->height;
    if (n>2)
    {
        for (int i = 0; i < n; i++)
        {
            matrix M(n-1);
            M = this->minor(0,i);
            d += pow(-1,i)*this->data[0][i] * M.determinant();
        }
        return(d);
    }
    else
    {
        d=(this->data[0][0])*(this->data[1][1])-(this->data[0][1])*(this->data[1][0]);
    }
    return(d);
}else{
    std::cout << "Matrix not squared!\n";
    return(d);
}
};

void matrix::determinant( double& determinant)
{
    determinant = this->determinant();
};

matrix matrix::inverse()
{//cout << "Inverse matrix started \t";
    matrix R(this->height);
    if (this->determinant() == 0)
    {std::cout << "Determinant = 0. No inverse matrix! Empty matrix returned! \n"; matrix O(this->height); return O;}
    else
    {
        for (int i = 0; i < this->height; i++)
            for (int j = 0; j < this->width; j++)
            {
                matrix M(height-1);
                M=this->minor(i,j);
                R.data[j][i] = pow(-1,i+j) * M.determinant();
            }
        R = R * (1/this->determinant());
    return R;
    }
};

void matrix::inverse_this()
{
    *this = this->inverse();
}

matrix matrix::elementwise_division(matrix B)
{
    if(height == B.height && width == B.width)
    {
        for(int i = 0; i < height; i++)
            {for(int j = 0; j < width; j++)
            {
                std::cout << data[i][j] << "\t" <<  B.data[i][j] << "\t\t";
                B.data[i][j] = data[i][j] / B.data[i][j];
                std::cout << B.data[i][j] << std::endl;
            }std::cout << std::endl;}
        return B;
    }else{
        std::cout << "matrix::elementwise_division(matrix B) Different sizes same matrix returned\n";
        return *this;
    }
};

int matrix::width_coord_of_max_elem_from_row(int column_index)
{
    int m_index = 0;
    /*if(this->width < column_index)
    {*/
    double m = 0;
    for(int i = 0; i < this->height; i++)
        if(fabs(this->data[i][column_index]) > m)
        {
            m = fabs(this->data[i][column_index]);
            m_index = i;
        }
    return m_index;
    /*}else{
    cout << "width_coord_max_elem_from_row left matrix ranges!\n";
    }
    return m_index;*/
};

void matrix::triangulation()
{//https://pro-prof.com/forums/topic/matrix-triangulation
    std::cout << "triangulation\n";
    if(this->height == this->width)
    {
        for(int i = 0 ; i < this->height; i++)
        {
            this->swap_rows(i, this->width_coord_of_max_elem_from_row(i));
            for(int j = i+1; j < this->height; j++)///присмотреться здесь
            {
                ///домножение строки на число и сложение строк
                for(int k = 0; k < this->height; k++)
                    this->data[k][j] -= this->data[k][i]*this->data[j][i]/this->data[i][i];
            }//cout << *this << endl;
            this->triangulation();
        }
    }else{
        std::cout << "Not square matrix!\n";
    }
}

void matrix::triangulation1(matrix &B)
{
const double elem_err = 0.0000000000001;
std::cout << "matrix triangulation with elem_err = " << elem_err << "\n";
if(width == height){

    if(width > 2)///при 2 рекурсивно не вызывается - обрывается
    {//cout << "matrix sizes>2x2!!!!\n";
        //cout << "in matrix\n" << *this << endl;
        //cout << "swap\n";
        ///находим строку с максимальным диагональным элементом и переставим столбцы
        int max_diag_elem = 0;
        for(int j = 0; j < 1; j++)
        {//cout << "j = " << j << endl;
            for(int i = j; i < height-1; i++)///самый большой диагональный элемент в строке наверх
            {//cout << "i = " << i << "\n";
                if(data[max_diag_elem][j] < data[i+1][j])
                    {max_diag_elem = i+1;std::cout << "new main element \t coordinate = " << max_diag_elem << "\telement = " << data[max_diag_elem][max_diag_elem] << std::endl;}
                //cout << data[max_diag_elem][max_diag_elem] << "\t" <<  data[i+1][i+1] << "\tmax el num =" << max_diag_elem << endl;
            }
            this->swap_rows(j, max_diag_elem);
            max_diag_elem = 0;
        }//cout << "swapped\n" << *this << endl;

        ///calc
        double m = 0;
        for(int i = 0+1; i < height; i++)///from second row
        {std::cout << "i = " << i << std::endl;
             m = -data[i][0]/data[0][0];
            for(int j = 0; j < width; j++)///for row
            {std::cout << "j = " << j << "\t";
                //if(data[0][0] != 0)
                //{

                    std::cout << "m = " << m << std::endl;
                    data[i][j] += data[0][j]*m;///minus norm first row
                    //B.data[j][0] = B.data[j][0]*(double)m;

                    if(fabs(data[i][j]) < elem_err)
                        data[i][j] = 0;
                //}
            }
        }
        matrix M(width-1);

        std::cout << "before triang\n" << *this << std::endl;
        M = this->minor(0,0);

        M.triangulation1(B);
        std::cout << "after triang\n" << M <<std:: endl;

        for(int i = 1; i < width; i++)
            for(int j = 1; j < width; j++)
                {
                    data[i][j] = M.data[i-1][j-1];
                    //cout << "index = " << i << "\t" << j << "\tdata = " << M.data[i-1][j-1] << endl;
                }
    }else{std::cout << "2x2!!!!\n";
    std::cout << *this << std::endl;
        double m = 0;
        m = -data[1][0]/data[0][0];
        std::cout << "m2 = " << m << std::endl;
        data[1][0] += data[0][0]*m;
        data[1][1] += data[0][1]*m;
    }
}else{std::cout << "Not squared matrix original returned\n";}
}

matrix matrix::average_column()
{
    double sum = 0;
    matrix r(1, width);
    for(int j = 0; j < width; j++)
    {
        for(int i = 0; i < height; i++)
            sum += data[i][j];
        r.data[0][j] = sum/height;
    sum = 0;
    }
    return r;
}

void matrix::column_multiply(int index_column, double m)
{
    if(index_column <= this->get_width())
    {
        for(int i = 0; i < this->get_height(); i++)
            this->set_data(i, index_column, this->get_data(i, index_column) * m);
    }else{
        std::cout << "matrix::column_multiply(int index_column, double m) -> matrix_width < index_column!, max index_column taken\n";
        this->column_multiply(this->get_width()-1, m);
    }
}

void matrix::column_shift(int index_column, double sh)
{std::cout << "column_shift()\n";
    if(index_column <= this->get_width())
    {
        for(int i = 0; i < this->get_height(); i++)
            this->set_data(i, index_column, this->get_data(i, index_column) + sh);
    }else{
        std::cout << "matrix::column_shift(int index_column, double sh) -> matrix_width < index_column!, max index_column taken\n";
        this->column_shift(this->get_width()-1, sh);
    }
}

matrix matrix::m_column_shift(int index_column, double sh)
{
    matrix r = *this;
        r.column_shift(index_column, sh);
    return r;
}

matrix matrix::matrix_exponential(int index_order)
{
    if(this->width == this->height)
    {
        matrix res, m_powered;
        double k_factorial = 1.0;

        ///index_order  == 0
        res.make_matrix_data(height, width, 1.0);

        if(index_order >= 1)
        {
            res += *this;
            k_factorial = 1.0;
            m_powered = *this;
        }

        for(int i = 2; i < index_order; i++)
        {
            m_powered *= *this;
            k_factorial *= i;
            res += m_powered/k_factorial;
        }
        return res;
    }
   std::cout << "matrix not squared 0 returned\n";
    return 0;
};

void matrix::Fibonachi2x2(int f_number)
{
    matrix base(2,2);
        base.set_data(0,0, 1.0); base.set_data(0,1, 1.0);
        base.set_data(1,0, 1.0); base.set_data(1,1, 0.0);
            std::cout << base << std::endl;

    if(f_number == 0)
        *this = base;

    matrix F(2,2); F = base;
    for(int i = 1; i < f_number; i++)
        {F *= base;std::cout << F << std::endl << std::endl;}
    *this = F;
};

///structure
void matrix::erase_row(int h)// - //minor 0,w
{
    if(h < height)
    {
    matrix result(height-1, width);
    for(int i = 0; i < height-1; i++)
    if( i < h )
    {
        for(int j = 0; j < width; j++)
            result.data[i][j] = data[i][j];
    }else{
        for(int j = 0; j < width; j++)
            result.data[i][j] = data[i+1][j];
    }
    *this = result;
    }else{
    std::cout << "matrix erase_row out of range matrix not changed!\n";
    }
}

void matrix::erase_column(int w)
{
    if(w < width)
    {
    matrix result(height, width-1);
    for(int i = 0; i < height; i++)
        for(int j = 0; j < width-1; j++)
            if(j < w)
            {
                result.data[i][j] = data[i][j];
            }else{
                result.data[i][j] = data[i][j+1];
            }
    *this = result;
    }else{
    std::cout << "matrix erase_column out of range matrix not changed!\n";
    }
};

matrix matrix::get_row(int h)
{
    matrix row(1, width);
    if( h < height)
    {
        for(int i = 0; i < width; i++)
            row.data[0][i] = data[h][i];
    }else{std::cout << "get_row out of range! zero returned\n";}
    return row;
}

matrix matrix::get_column(int w)
{
    matrix column(height, 1);
    if( w < width)
    {
        for(int i = 0; i < height; i++)
            column.data[i][0] = data[i][w];
    }else{std::cout << "get_column out of range! zero returned\n";}
    return column;
}

matrix matrix::get_matrix_part(int h1, int h2, int w1, int w2)///from h1,w1 to h2,w2
{//cout << "get_matrix_part(int h1, int h2, int w1, int w2)\t" << h1 << " " << h2 << " " << w1 << " " << w2 << endl;
    //this->info();
    ///prep
    if(h1!=0 || h2!=0 || w1!=0 || w2!=0)///если процедура не требуется
    {
        if(height != 0)
        {
            int k;
            if(h1 <= 0) h1 = 0;
            if(h2 <= 0) h2 = 0;
            if(h1 > h2)
                {k = h2;h2 = h1;h1 = k;}
            if(h2 > height) h2 = height;

            if(w1 <= 0) w1 = 0;
            if(w2 <= 0) w2 = 0;
            if(w1 > w2)
                {k = w2;w2 = w1;w1 = k;}
            if(w2 > width) w2 = width;
        //cout << h1 << " " << h2 << " " << w1 << " " << w2 << endl;
            //if(w1 > width-1) w1 = width-1;
        ///func
        matrix r(h2-h1, w2-w1); //cout << h1 << "\t" << h2 << "\t" << w1 << "\t" << w1 << endl;
        for(int i = h1; i < h2; i++)
            for(int j = w1; j < w2; j++)
               r.set_data(i - h1, j - w1, data[i][j]);
        return r;
        }
    }
    return *this;
};

matrix matrix::merge_height(const matrix B, int pozition_row)
{
    if(height == 0)
        return B;
    if(B.height == 0)
        return *this;

    if(width == B.width)
    {
        if(pozition_row > height) /// равносильно добавлению в конец в следующем if
                pozition_row = height;

        if(pozition_row >= 0)
        {
            matrix AB(height + B.height, width);

            for(int j = 0; j < width; j++)
                for(int i = 0; i < pozition_row; i++)
                    AB.data[i][j] = data[i][j];

            for(int j = 0; j < width; j++)
                for(int i = pozition_row; i < B.height + pozition_row; i++)
                    AB.data[i][j] = B.data[i-pozition_row][j];

            for(int j = 0; j < width; j++)
                    for(int i = B.height + pozition_row; i < height + B.height; i++)
                        AB.data[i][j] = data[i-B.height][j];

            return AB;
        }else{
        std::cout << "merge_height pozition_row < 0! Same matrix returned\n";
        return *this;
        }
    }else{
        std::cout << "merge_height different width of matrixes! Same matrix returned\n";
        return *this;
    }
};

matrix matrix::merge_height(matrix B)///end
    {return this->merge_height(B, this->get_height());}

matrix matrix::merge_width(matrix B, int pozition_column)
{
    if(height == 0)
        return B;
    if(B.height == 0)
        return *this;

    if(height == B.height)
    {
        if(pozition_column > width) /// равносильно добавлению в конец в следующем if
                pozition_column = width;

        if(pozition_column >= 0)
        {
            matrix AB(height, width + B.width);

            for(int i = 0; i < height; i++)
            {
                for(int j = 0; j < pozition_column; j++)
                    AB.data[i][j] = data[i][j];

                for(int j = pozition_column; j < pozition_column + B.width; j++)
                    AB.data[i][j] = B.data[i][j-pozition_column];

                for(int j = pozition_column + B.width; j < B.width + width; j++)
                    AB.data[i][j] = data[i][j-B.width];
            }
            return AB;
        }else{
        std::cout << "merge_width pozition_column < 0! Same matrix returned\n";
        return *this;
        }
    }else{
        std::cout << "merge_width different height of matrixes! Same matrix returned\n";
        return *this;
    }
};

matrix matrix::merge_width(matrix B)///end
    {return this->merge_width(B, this->get_width());}

void matrix::merge_height(matrix A, matrix B, int pozition_row)
{//cout << "MERGE!!!\n" ;//<< pozition_row << endl;
    if(A.height == 0 && A.width == B.width)
    {//cout << "row pozition = " << pozition_row << endl;
        *this = B;
    }
    if(B.height == 0 && A.width == B.width)
    {//cout << "row pozition = " << pozition_row << endl;
        *this = A;
    }
    if(A.height == 0 && A.width == 0)
    {//cout << "row pozition = " << pozition_row << endl;
        *this = B;
    }
    ///A.width - >A.height
    if(pozition_row <=A.height && pozition_row >= 0 && A.height != 0 && B.height != 0)
    {//cout << "row pozition = " << pozition_row << endl;
    matrix AB(A.height + B.height, A.width);

    for(int j = 0; j < A.width; j++)
        for(int i = 0; i < pozition_row; i++)
            AB.data[i][j] = A.data[i][j];

    for(int j = 0; j < A.width; j++)
        for(int i = pozition_row; i < B.height + pozition_row; i++)
            AB.data[i][j] = B.data[i-pozition_row][j];

    for(int j = 0; j < A.width; j++)
        for(int i = B.height + pozition_row; i < A.height + B.height; i++)
            AB.data[i][j] = A.data[i-B.height][j];

    *this = AB;
    }///else{cout << "Merge_height wrong position_column!\n";}
};


void matrix::merge_width(matrix A, matrix B, int pozition_column)
{
    if(pozition_column <=A.width && pozition_column >= 0)
    {
    matrix AB(A.height, A.width+B.width);
    for(int i = 0; i < A.height; i++)
    {//cout << "fori = " << i << endl;
        for(int j = 0; j < pozition_column; j++)
            AB.data[i][j] = A.data[i][j];
        //cout << "0\n" << AB << endl;
        for(int j = pozition_column; j < pozition_column + B.width; j++)
            AB.data[i][j] = B.data[i][j-pozition_column];
        //cout  << "pc\n" << AB << endl;
        for(int j = pozition_column + B.width; j < B.width + A.width; j++)
            AB.data[i][j] = A.data[i][j-B.width];
        //cout << "pc + Bw\n" << AB << endl;
    }
    *this = AB;
    }else{std::cout << "Merge_width wrong position_column!\n";}
};

void matrix::swap_elements(int i1, int j1, int i2, int j2)
{
    if( i1 < this->height && i2 < this->height && j1 < this->width && j2 < this->width)
    {
        double buff = 0;
        buff = this->data[i1][j1];
        this->data[i1][j1] = this->data[i2][j2];
        this->data[i2][j2] = buff;
    }else{
    std::cout << "swap_elements left matrix ranges - original returned\n";
    }
}

void matrix::swap_columns(int w1, int w2)
{
    if( w1 < this->width && w2 < this->width)
    {
    for(int i = 0; i < this->height; i++)
        swap_elements(i, w1, i, w2);
    }else{
    std::cout << "swap_columns left matrix ranges - original returned\n";
    }
};

void matrix::swap_rows(int h1, int h2)
{
    if( h1 < this->height && h2 < this->height )
    {
        if(h1 != h2)
            for(int i = 0; i < this->width; i++)
                swap_elements(h1, i, h2, i);
    }else{
    std::cout << "swap_rows left matrix ranges - original returned\n";
    }
};

///types

/*matrix matrix::array_to_matrix_column(double* arr, int length)
{
    matrix r(length, 1);
    for(int i = 0; i < r.get_height(); i++)
        r.set_data(i, 0, arr[i]);
    return r;
};

matrix matrix::array_to_matrix_row(double* arr, int length)
{
    matrix r(1, length);
    for(int i = 0; i < r.width(); i++)
        r.set_data(0, i, arr[i]);
    return r;
};*/

void matrix::array_to_matrix_column(double* arr, int length)
{
    matrix r(length, 1);
    for(int i = 0; i < r.get_height(); i++)
        r.set_data(i, 0, arr[i]);
    *this = r;
};

matrix matrix::matrix_column_from_array(double* arr, int length)
{
    matrix r(length, 1);
    for(int i = 0; i < r.get_height(); i++)
        r.set_data(i, 0, arr[i]);
    return r;
};

void matrix::array_to_matrix_row(double* arr, int length)
{
    matrix r(1, length);
    for(int i = 0; i < r.get_width(); i++)
        r.set_data(0, i, arr[i]);
    *this = r;
};

matrix matrix::matrix_row_from_array(double* arr, int length)
{
    matrix r(1, length);
    for(int i = 0; i < r.get_width(); i++)
        r.set_data(0, i, arr[i]);
    return r;
};

///FILES
void matrix::Find_Size_of_matrix_file(const std::string file_name)
{
    std::string buff = "";
    std::ifstream fin(file_name);
    height = 0; width = 0;
    while(!fin.eof())
    {
            getline(fin, buff);
            height++;
    }
    fin.close();//cout << "height = " << height << "\t";

    fin.open(file_name);
    while(!fin.eof())
    {
            fin >> buff;
            width++;
    }
    width = (int)width/height;
    fin.close();//  cout << "width = " << width << endl;
}

void matrix::Find_Size_of_matrix_emptstringend_file(const std::string file_name, int emptstringend_num)
{
    std::string buff = "";
    std::ifstream fin1(file_name);
    height = 0; width = 0;
    while(!fin1.eof())
    {
            getline(fin1, buff);
            height++;
    }
    height-=emptstringend_num;
    fin1.close();//cout << "height = " << height << "\t";

    std::ifstream fin2(file_name);
    while(!fin2.eof())
    {
        fin2 >> buff;
            width++;
    }//cout << "elcount = " << width << endl;
    width = width/height;
    fin2.close();//cout << "width = " << width << endl;
}

void matrix::Load_matrix(const std::string file_name)
{std::cout << "Load matrix started\t" << file_name <<"\t";
    if(FileExists(file_name.c_str()))
    {
        Find_Size_of_matrix_file(file_name);//cout << "1!!!!!!!!!!!!\n";
        std::ifstream fin(file_name);
        data = new double* [height];
        for(int i = 0; i < height; i++)
        {
            data[i] = new double[width];
            for(int j = 0; j < width; j++)
                fin >> data[i][j];
        }
        fin.close();
    }else{
        matrix m_empty;
        *this = m_empty;
        std::cout << "WARNING File do not exist! \t" << file_name << "\tEmpty matrix returned\n";
    }
    std::cout << "Matrix loaded\n";
}
void matrix::Load_matrix_emptstringend(const std::string file_name, int emptstringend_num)
{std::cout << "Load matrix started\t" << file_name <<"\t";
    if(FileExists(file_name.c_str()))
    {
        Find_Size_of_matrix_emptstringend_file(file_name, emptstringend_num);//cout << "1!!!!!!!!!!!!\n";
        std::ifstream fin(file_name);
        data = new double* [height];
        for(int i = 0; i < height; i++)
        {
            data[i] = new double[width];
            for(int j = 0; j < width; j++)
                fin >> data[i][j];
        }
        fin.close();
    }else{
        matrix m_empty;
        *this = m_empty;
        std::cout << "No such file! \t" << file_name << "\tEmpty matrix returned\n";
    }
    std::cout << "Matrix loaded\n";
}

void matrix::Save_matrix(const std::string file_name)
{std::cout << "Save matrix started\t" << file_name <<"\t";
    std::ofstream fout(file_name);
        fout << *this;
    fout.close();
    std::cout << "Matrix saved\n";
}
void matrix::Save_matrix_end(const std::string file_name)
{
    std::cout << "Save matrix end started\t" << file_name <<"\t";
    std::ofstream fout(file_name, std::ios::app);
        fout << *this << std::endl;
    fout.close();
    std::cout << "Matrix saved\n";
};


void matrix::info()
{
    std::cout << "Matrix object\t height = " << this->height << "\twidth = " << this->width << std::endl;
};

bool FileExists(std::string fname)///old const char *fname - > string
{
    //return ::ifstream(fname, ios::in ) != NULL;////::ifstream(fname, ios::in | ios::nocreate) != NULL;
    std::ifstream f(fname.c_str());
    return f.good();
}

///io.h

bool FileExists1(std::string fname)
{
    std::ifstream file;
    file.open(fname);
    file.close();
    if (file)
    {
        std::cout << "Success.\n";
        return true;
       /*int n;
        file >> n; // oops
        std::cout << n;*/
    }
    return false;
}
