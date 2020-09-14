#ifndef DIFEQ_CALC_H_INCLUDED
#define DIFEQ_CALC_H_INCLUDED

///Equations
void Calc_VdP_equation_4(std::string Data_Way, std::string base_file_name, int fix_max);
void Calc_VdP4_FrontRear(std::string Data_Way, std::string base_file_name, std::string mode, int fix_max);///30.12.19
void Calc_VdP_equation_4_mdyn(std::string Data_Way, std::string base_file_name, std::string mode, int fix_max);

#endif // DIFEQ_CALC_H_INCLUDED
