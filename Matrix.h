#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <iostream>


class matrix
{
private:
protected:
    int height;
    int width;
    double** data; ///можно перевести основу на vector
public:

///constructors//des//copy
    matrix(); ///empty
    void matrix_alloc(int , int); ///создание пустых массивов размера
    matrix(int , int); ///rectangle
    matrix(int); ///Square matrix
    matrix(int n, char);///квадратная 'E' - единичная, 'n' - в столбцах натуральные числа с 1 по n.
    matrix(int i, int j, char);///'E' - единичная, 'n' - в столбцах натуральные числа с 1 по n.
    ~matrix(); ///destructor
    matrix(const matrix&); ///copy

///primary
    void delete_data();///удаляет данные матрицы заменяя пустой с размером 0х0 /// можно оптимизировать
    int get_height() const;
    int get_width() const;
    double get_data(int i, int j) const;
    void set_height(int h);
    void set_width(int w);
    void set_data(int i, int j, double d);

    void set_column_data(int h, int w, double data);///заполняет указанный столбец с 0 до h, стоблец с индексом w, если нет стобца - создает.

    void make_matrix_data(int h, int w, double d);

    void insert_array_to_matrix_column(int index_column, double* arr, int length);///overwrite and fit to matrix.height();
    void insert_array_to_matrix_row(int index_row, double* arr, int length);///overwrite and fit to matrix.width();

    double get_max_el_column(int index_column);
        double get_max_el_column_gap(int index_column, int gap_index_from, int gap_index_to);
    double get_min_el_column(int index_row);

    int get_max_elind_column_gap(int index_column, int gap_index_from, int gap_index_to);
    int get_min_elind_column_gap(int index_column, int gap_index_from, int gap_index_to);

///operatros
    matrix& operator=(const matrix& A);

    matrix& operator+=(const matrix& A);
    friend matrix operator+(const matrix& A, const matrix& B);

    matrix operator* (double a); /// можно оптимизировать /// реализвоать ч/з Алгоритм Штрассена O(n^2.81)
    friend matrix operator*(const double& a, const matrix& A);
    friend matrix operator*(const matrix& A, const matrix& B);
    matrix operator*=(const matrix& A);

    matrix operator/ (const double a);

    matrix operator-= (const matrix& A);
    friend matrix operator-(const matrix& A, const matrix& B);

    friend std::ostream& operator << (std::ostream &C, const matrix &A);
    friend std::istream& operator >> (std::istream &C, matrix &A);

///FUNCTIONS
///math
    matrix abs_column(int index_column);
    void abs_index_column(int index_column);
    //matrix abs_row( int index_row);
    //matrix abs();///all

    matrix transpose();
        void transpose(matrix& ); /// как поступать в таких случаях?
    matrix minor(int l,int k); ///"вычеркнуть" номер столбца и строки /// квадратная матрица!!!!!
    double determinant(); ///default algorithm  /// minor /// можно оптимизировать
        void determinant(double& determinant);
    matrix inverse();
        void inverse_this();
    matrix elementwise_division(matrix B);/// sizes1==sizes2 /// A/B == A.elementwise_division(B)

    int width_coord_of_max_elem_from_row(int column_index);///coordinates counts from 0 //?
    void triangulation(); ///?
    void triangulation1(matrix& B); ///elem err ///?

    matrix average_column();

    void column_multiply(int index_column, double multiplier);
    void column_shift(int index_column, double sh);///прибавить к каждому элемента столбца index_column sh
        matrix m_column_shift(int index_column, double sh);///производная функция

    matrix matrix_exponential(int index_order);

///structure
    void erase_row(int h);/// - /// удаляет строку в данной матрице /// можно оптимизировать
    void erase_column(int w);///| /// удаляет столбец в данной матрице /// можно оптимизировать
    matrix get_row(int h);/// доставляет строку данной матрицы по номеру /// можно оптимизировать
    matrix get_column(int w);
    matrix get_matrix_part(int h1, int h2, int w1, int w2); ///simple part from h1,w1 to h2,w2 ///доставляет часть матрицы начиная со столбца h1 до h2 и со строки w1 до w2

    matrix merge_height(const matrix B, int pozition_row); ///отлажено
        matrix merge_height(matrix B);///в конец
    matrix merge_width(matrix B, int pozition_column); ///отлажено
        matrix merge_width(matrix B);///в конец
        void merge_height(matrix A, matrix B, int pozition_row); /// НЕ РАБОТАЕТ!!!
        void merge_width(matrix A, matrix B, int pozition_column); //

    void swap_elements(int i1, int j1, int i2, int j2);///можно сделать через указатели /// можно оптимизировать
        void swap_columns(int w1, int w2);
        void swap_rows(int h1, int h2);

///types
    //matrix array_to_matrix_column(double* arr, int length);
   // matrix array_to_matrix_row(double* arr, int length);
    void array_to_matrix_column(double* arr, int length);
    matrix matrix_column_from_array(double* arr, int length);
    void array_to_matrix_row(double* arr, int length);
    matrix matrix_row_from_array(double* arr, int length);

///files
    void Find_Size_of_matrix_file(const std::string file_name); ///устанавливает height и width данной матрицы по размерам строк и столбца, конец файла на числе /// Можно оптимизировать
    void Find_Size_of_matrix_emptstringend_file(const std::string file_name, int emptstringend_num);///когда в конце файла пустая строка
    void Load_matrix(const std::string file_name); ///memory allocation
    void Load_matrix_emptstringend(const std::string file_name, int emptstringend_num);
        //void Load_matrix_with_header(const string file_name);///сделать /// вынести в отдельный блок
    void Save_matrix(const std::string file_name);
    void Save_matrix_end(const std::string file_name); /// сохранение матрицы в конец файла

///extra
    void info();
};
///outer
    bool FileExists(std::string fname);
    bool FileExists1(std::string fname);
    bool FileExists2(const char *fname);

#endif // MATRIX_H_INCLUDED
