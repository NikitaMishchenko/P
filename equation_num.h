#ifndef EQUATION_NUM_H_INCLUDED
#define EQUATION_NUM_H_INCLUDED

#include <functional>



class DifEq
{
private:
protected:
public:
    //virtual DifEq() = 0;
    //virtual ~DifEq() = 0;

    void calculate();
};

class DifEq_2order : public DifEq
{
private:
    ///coefficients for equation d2y/dx2 + f(x)dy/dx+g(x)=0
    double fa0,fa2,fa4; ///f(x) = fa0x^0+fa2x^4+fa4x^4
    double ga1,ga3,ga5; ///f(x) = ga0x^1+ga2x^3+ga4x^5
    matrix result; ///time, y(x), z(y), lam1Re, lam1Im, lam2Re, lam2Im

protected:
public:
    //DifEq_2order() override;
    //~DifEq_2order() override;

    void calculate();
};

void base_equation(std::string , std::string, double, double);

double p(double, double, double, double);
double p1(double, double, double, double, double);
double g(double, double, double, double);
double g1(double, double, double, double, double);

matrix find_incline_mz_points(std::string Data_Way, std::string base_file_name);

matrix equation_func1(std::string Data_Way, std::string file_name, double h, int N, double gb, double gc, double pb, double pc, double x0, double y0, double z0, double shift_angle);///система линейных уравнений
///система линейных уравнений
matrix equation_func2(std::string Data_Way, std::string file_name, double h, int N, double gb0, double gb2, double gc, double pb, double pc, double x0, double y0, double z0, double shift_angle);///система линейных уравнений
matrix equation_func3(std::string Data_Way, std::string file_name, double h, int N,int fix_max, double gb0, double gb1, double gb2, double gc, double pb, double pc, double x0, double y0, double z0, double shift_angle);///система линейных уравнений
matrix equation_func4(std::string Data_Way, std::string file_name, double h, int N,int fix_max, double gb0, double gb2, double gb4, double gc, double pb, double pc, double x0, double y0, double z0, double shift_angle);///система линейных уравнений
matrix equation_func_VdP_4(std::string Data_Way, std::string file_name, double h, int N,int fix_max, double gb0, double gb2, double gb4, double gc, double pb, double pc, double x0, double y0, double z0);
matrix equation_func_VdP_4_fr(std::string Data_Way, std::string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4, double gc, double pb, double pc, double fr, double x0, double y0, double z0);
matrix equation_func_VdP_4_mdyn(std::string Data_Way, std::string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4, double gc, double gd, double pb, double pc, double pd, double mdyn, double x0, double y0, double z0);
matrix equation_func_VdP_4_mdyn_poly_mst(std::string Data_Way, std::string file_name, double h, int N, int fix_max, double gb0, double gb2, double gb4, double gc0, double gc2, double gc4, double gd, double pb, double pc, double pd, double mdyn, double x0, double y0, double z0);
matrix equation_func_polynomial4order_f_g(std::string Data_Way, std::string file_name, double h, int N, int fix_max,
                                          std::function<double(double fa0, double fa2, double fa4, double f_arg)> f_poly,
                                          std::function<double(double ga0, double ga2, double ga4, double g_arg)> g_poly,
                                          double gb0, double gb2, double gb4,double gc0, double gc2, double gc4,
                                          double gd, double pb, double pc, double pd,double mdyn, double x0, double y0, double z0);

void eq_report();///нужно выносить в отдельный класс

#endif // EQUATION_NUM_H_INCLUDED
