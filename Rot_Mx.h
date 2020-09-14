#ifndef ROT_MX_H_INCLUDED
#define ROT_MX_H_INCLUDED

#include "string.h"

class Rot_Mx
{
private:
protected:

    matrix freq;

public:
    matrix voltage;///time_volt, volt
    //matrix time;
    matrix velocity;///time_v, v
    double time_shift; ///
    std::string holl_conf;///количество и конфигурация магнитов

    Rot_Mx();
    Rot_Mx(int h);

    Rot_Mx(std::string file_name);

    void set_holl_conf(std::string h_c);
    std::string get_holl_conf();

    //void Calc_velocity();///из данных напряжения датчика холла строиятся времена прохождений - скорость
    matrix Get_Calc_velocity();///из данных напряжения датчика холла строиятся времена прохождений - скорость ///settings file
    matrix Get_Calc_velocity_slope(int window_range); /// "approximation.h" -> Find_Incline_Coefficient(Part, &sl, window_range);

    matrix Get_Avg_velocity(Rot_Mx A, double from_v);///current + A from from_v till end


    void Load_Rot_Mx(std::string file_name);///построение по напряжению
    void Load_Rot_Mx_velocity(std::string file_name); /// построение по посчитаной частоте (velocity)
    void Save_Rot_Mx(std::string file_name);

    void info();
};
double get_time_shift(matrix Data, double f);

matrix Rot_Mx_time(std::string Data_Way, std::string data_file_name, double ADC_freq);
void Rot_Mx_test(std::string Data_Way, std::string data_file_name, std::string h_c, double ADC_freq, int slope_window_range, int slope2_window_range); ///h_c - конфигурация датчика холла

#endif // ROT_MX_H_INCLUDED
