#ifndef OSC_CLASS_H_INCLUDED
#define OSC_CLASS_H_INCLUDED

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "matrix.h"
#include "Approximation.h"

class oscillation ///counter / 3flow /angle /dangle /time
{///?лучше ли реализовать как вектор структуры (int, flow (double, double, double), angle, dangle, time)?
private:
    double init_angle;
protected:
    matrix time; ///явно в angle и flow
    matrix angle; ///зарегистрированный угол атаки от вреени
    double calibr_angle;
    matrix dangle; /// производная угла по времени
    matrix flow; /// потоковые данные
    matrix freq; /// частота колебаний
    matrix balance_angle;

    std::string model_name;///название модели
    std::string data_source;
    double I;///момент инерции
    double S;///характерная площадь
    double l;///характерная длина

public:
    oscillation(); ///empty
    oscillation(int h);
    oscillation(std::string file_name); /// приведенный по формату файл

    ///main
    //matrix count_freq(); ///3p-derevative
    double calc_balance_angle();

    matrix count_freq(matrix t, matrix A);
    void callibrate_to_zero();
    void callibr_angle();

    void shift_time(double time_shift);
    void shift_angle(double angle_shift);
    void angle_shift_balance(int);
    void angle_shift_balance();

    void callibrate_to_radian();
    void callibrate_to_degree();

    ///primary
    double get_init_angle();
    //void set_balance_angle(double);
    //double get_balance_angle();

    void set_time(matrix );
    matrix get_time();

    void set_angle(matrix );
    void set_angle_data(int, int, double);
    matrix get_angle();

    void set_calibr_angle(double a);
    double get_calibr_angle();

    void set_dangle(matrix );
    matrix get_dangle();

    void set_freq(matrix f);
    matrix get_freq();

    void set_flow(matrix );
    matrix get_flow();
    double get_flow_Re(int i);///Reynolds num in 10^6
    double get_flow_M(int i);///Mach num
    double get_flow_v(int i);///veloc
    double get_flow_q(int i);///head in Pa


    void set_flow_data(int i, int j, double data);
    double get_flow_data(int i, int j);

    void set_balance_angle(matrix );
    void set_balance_angle_data(int , int, double);
    void set_balance_angle_data(int , double);
    matrix get_balance_angle();
    double get_balance_angle_data(int i, int j);
    double get_balance_angle_data(int i);

    void set_data_source(std::string );
    std::string get_data_source();

    void set_model(std::string );
    std::string get_model();
    void set_model_I(double );
    double get_model_I();
    void set_model_S(double );
    double get_model_S();
    void set_model_l(double );
    double get_model_l();

    ///New Structures
    void cut_from_to(double time1, double time2);
    oscillation oscillation_cut_from_to(double time1, double time2);


    matrix w_harm();
    matrix w_harm1();
    void middle_angle_timeavg_halfT(const matrix& env_top, const matrix& env_bot, matrix& middle_angle);
    matrix middle_angle_calc(matrix env_top, matrix env_bot);

    ///stream
    friend std::ostream& operator << (std::ostream &C, const oscillation &A);
    ///oscillation& operator=(const oscillation& A);


    ///Data save/load
    ///SAVE
    void Load_Oscillation_Raw_Protocol(const std::string Data_Way, const std::string base_file_name);
    double Load_get_Starting_Impolse_Time(std::string Data_Way, std::string base_file_name);
    matrix Flow_from_PTL(std::string Data_Way, std::string base_file_name); /// g = 9.80665;
    int Load_get_Starting_Impolse_Number(std::string Data_Way, std::string base_file_name);
    void Load_Oscillation_flow_Protocol(std::string Data_Way, std::string base_file_name, double t_shift);
    void Load_Oscillation_caliber_angle(std::string Data_Way, std::string base_file_name);
    void Load_Oscillation_Protocol(std::string Data_Way, std::string base_file_name);
    void Load_Oscillation_Prop_File(std::string Data_Way, std::string base_file_name, std::string mode);
    void Load_Oscillation_Prop_File(std::string Data_Way, std::string base_file_name);
    void Load_flow_summary(std::string DataWay, std::string file_name);///вставляет название модели если не было до этого
    void Load_model_info(std::string DataWay, std::string file_name);
    void Load_balance_angle_summary(std::string DataWay, std::string file_name);///по имени протокола = data_source через временные окна
    ///SAVE
    void Save_Oscillation(std::string Data_Way, std::string base_file_name);

    ///
    void info();
};

#endif // OSC_CLASS_H_INCLUDED
