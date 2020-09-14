#ifndef MANAGER_H_INCLUDED
#define MANAGER_H_INCLUDED

#include <iostream>
#include <queue>
///#include <functional>

#include "matrix.h"
#include "Osc_Class.h"
#include "osc.h"

class manager
{
private:
    bool sucess = false;
    const std::string mark_file_manager_begin = "MANAGER";
    const std::string mark_file_manager_end = "END_MANAGER";
    const std::string mark_file_commands = "commands:";
    const std::string mark_file_proceded_commands = "proceded_commands:";
    const std::string mark_file_data_paths = "data_paths:";
    const std::string manager_file_format = ".txt";
    const std::string own_data_path;
protected:
public:
    std::vector <std::string> commands; ///набор команд к исполнению
    std::vector <std::string> proceded_commands; /// набор выполненных комманд /// в private
    std::vector <std::string> data_paths; /// пути к данным
    std::string stop_command;
    std::string name;

    manager();
    ~manager();

    void manage_commands(); /// выполнение всех команд
    void manage_command_index(int index_command); /// выполнение конкретной комманды
    void is_in_commands(std::string command); ///есть ли в списке комманд
    void is_in_proceded_commands(std::string pcommand); ///есть ли в списке комманд
    bool check_proceded_commands();///все ли команды выполнены?
    void reproceed();///выполнить команды до соответстви€ заданным, если это возможно
    bool reproceedable();

///GET/SET/ADD/ERASE
    ///comamnds
    void set_commands(std::vector <std::string> new_commands);
    void push_command(std::string new_command);
    ///proceded_commands
    void set_proceded_commands(std::vector <std::string> new_p_commands);
    void push_proceded_command(std::string new_p_command);


///FILES
    void Load_manager(std::string file_name);
    void Save_manager(std::string file_name);

///EXTRA
    std::string get_info();
        void info();
    void info_load_from_file();
    void create_load_from_file_fishfile(const std::string full_file_name);
};

std::string manager_int_num(int i);

std::string to_string(double r);
std::string my_to_string(bool b);

///GNUPLOT assosiation


///OSILLATION
struct manager_instruction
{
    std::string summary_file;
    std::string DataWay;
    std::string base_file_name;
};
std::queue<manager_instruction> manager_read(std::string full_manager_file_name);
void manager_amplitude_analysis(std::queue<manager_instruction>);
void General_Procedure_Oscillation(std::string DataWay_summary, std::string Data_Way, std::string base_file_name, std::string name_add, std::string mode_middle_angle_correction);


///Fluent report
void procedure_get_my_res(std::string Data_Way, std::string Data_Way_Res, std::string base_file_name, double angle);
void Flow_pulsation_check(std::string Data_Way, std::string file_name);

///INERTIA
void Fix_LIR_Data(std::string Data_Way, std::string base_file_name, std::string result_file, double border); ///исправление данных, при зашкаливании данных угла (например, вращение)

void procedure_inertia(std::string Data_Way, std::string base_file_name, int  N_p, int overlap, int fft_zero_index, double sigma, int K_pow, double m1, double r1, double g); ///гаусс-окно
double cone_inert_2810(double h, double m2, double xi, double alpha1, double T);
double inertia_cg(double xl, double yl, double h, double m1, double m2, double alpha1, double alpha3);
double inertia_cg1(double xl, double yl, double h, double m1, double m2, double a_free, double a_aft, double a_w);

double inertia_momentum(double xl, double yl, double h, double m1, double m2, double a_free, double a_aft, double a_w, double g, double T);///«ƒ≈—№ Ќ≈ ”„»“џ¬ј≈“—я m1 !!!
void procedure_inertia_analysis(double Isw, double m1, double m2, double h, double d, double alpha0, double alpha1, double alpha3, double yl, double r2, double g, double T);

matrix procedure_get_pendulum_frequency(std::string Data_Way, std::string base_file_name, oscillation D, std::string method_find_freq, int window_size, int step_size, int zero_index, int K_power, double avg_angle);///time // freq ///f = 1/T
matrix procedure_get_pendulum_frequency(std::string Data_Way, std::string base_file_name, oscillation D, std::string method_find_freq, int window_size, int step_size, int zero_index, int K_power);///time // freq ///f = 1/T
void procedure_moment_of_inertia_stand(std::string Data_Way, std::string base_file_name, std::string method_find_freq, int window_size, int step_size,int zero_index, int K_power);
matrix procedure_Redd(matrix envelop, oscillation D);
void procedure_ampl_analisys(std::string Data_Way, std::string base_file_name, std::string name_add, oscillation& D, double freq, double angle_shift);
void procedure_ampl_analisys(std::string Data_Way, std::string base_file_name, oscillation& D, double freq, double angle_shift);
matrix data_cutter(matrix A, int column_num, double column_data1, double column_data2);
void Redd_Cut_file(std::string Data_Way, std::string instruct_name);

matrix procedure_osc_analisys(std::string Data_Way, std::string base_file_name, oscillation D, int window_size, int step_size, int zero_index);
///FFT
void procedure_Mz_FFT(std::string Data_Way, std::string base_name, double time_from, double time_to, double window_step_time, int window_size_power, int window_zeroes_power);

///Equations
//void Calc_VdP_equation_4(string Data_Way, string base_file_name, int fix_max);
//void Calc_VdP4_FrontRear(string Data_Way, string base_file_name, string mode, int fix_max);///30.12.19
//void Calc_VdP_equation_4_mdyn(string Data_Way, string base_file_name, string mode, int fix_max);

void Analys_VdP_equation_4(std::string Data_Way, std::string base_file_name, std::string name_add,double mu);
void Analys_VdP4_mdyn(std::string Data_Way, std::string base_file_name, std::string name_add);
void Analys_VdP4_mdyn(std::string Data_Way, std::string base_file_name, std::string name_add, std::string angle_val, double mdyn, double epsilon4, double epsilon2 );
void Calc_Equation_procedure();
void Get_Cycle_Radii(std::string Data_Way, std::string file_name, double from, double to, double step);


///FLUENT REPORT
void get_mz_results(std::string Data_Way, std::string base_name, std::string result_name, int from_index, int to_index);
void get_mz_results1(std::string Data_Way, std::string base_name, std::string result_name, int from_index, int to_index);
void symmetry_coeff_report_flu(std::string Data_Way, std::string file_name);
#endif // MANAGER_H_INCLUDED
