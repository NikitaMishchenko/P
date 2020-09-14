#include "manager.h"

#include <iostream>
#include <cmath>

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <string>
#include <queue>
#include <thread>
#include <mutex>

#include <functional>


#include "matrix.h"
#include "Approximation.h"
#include "Osc_Class.h"

#include "osc.h"
#include "Rot_Mx.h"
#include "FFT.h"
#include "equation_num.h"

using namespace std;

const double PI = 3.1415926535897932384626433832795;


manager::manager()
{
    sucess = false;
    stop_command = "STOP";
    name = "";
    vector <string> n_c;
    vector <string> n_pc;
    vector <string> n_dp;
        commands = n_c;
        proceded_commands = n_pc;
        data_paths = n_dp;
}

manager::~manager()
{
    sucess = false;
    commands.clear();
    proceded_commands.clear();
    data_paths.clear();
}

void manager::manage_commands()
{cout << "manager::manage_commands()" << endl;
    cout << "\t";
    for (auto c : commands)
    {
        ///выполнение комманд
        cout << c << "\t";

        if(true && c != stop_command)///если не было фэйла при выполнении команды
            sucess = true;
        ///запись выполненных команд
        proceded_commands.push_back(c);
        if(c == stop_command)
            break;
    }
    cout << endl;
}

void manager::manage_command_index(int index_command)
{cout << "manager::manage_commands()" << endl;
    string command;
        command = commands.at(index_command);
    proceded_commands.push_back(command);
}

bool manager::check_proceded_commands()
{cout << "manager::check_proceded_commands() = ";
    if( commands == proceded_commands)
    {
        return true;
        cout << endl;
    }cout << endl;
    return false;
};

bool manager::reproceedable()
{cout << "bool manager::reproceedable()" << endl;
    if(proceded_commands.size() == 0)///можно выполнять любые команды
        return false;
    if(commands.size() >= proceded_commands.size())
    {//cout << "\tif\t" << commands.size() << "\t" << proceded_commands.size() << endl;
        vector <string>::const_iterator it_c = commands.begin(), it_pc = proceded_commands.begin();
        while(it_pc != proceded_commands.end())/// если все сделанные команды совпадают, то все ок
        {//cout << "\t\twhile: \t" << *it_c << " ?= " << *it_pc << endl;
            if(*it_c == *it_pc)
            {
                it_c++;it_pc++;
            }else{
                return false;
            }
        }//cout << "REPROCEEDEBLE\n";/// можно выполнять несделанные команды
        return true;
    }else{///выполненных команд больше чем заданных
        cout << "WARNING commands.size() > proceded_commands()" << commands.size() << "\t" << proceded_commands.size() << endl;
        return false;
    }
    cout << "WARNING smth happen false returned\n";
    return false;
}

///GET/SET/ADD/ERASE
void manager::set_commands(vector <string> new_commands)
{
    commands = new_commands;
    sucess = this->check_proceded_commands();
};

void manager::push_command(string new_command)
{
    commands.push_back(new_command);sucess = false;
    sucess = this->check_proceded_commands();
};

void manager::set_proceded_commands(vector <string> new_p_commands)
{
    proceded_commands = new_p_commands;
    sucess = this->check_proceded_commands();
};
void manager::push_proceded_command(string new_p_command)
{
    proceded_commands.push_back(new_p_command);
    sucess = this->check_proceded_commands();
};

///FILES
void manager::Load_manager(string full_file_name)
{cout << "void manager::Load_manager(string full_file_name)" << endl;
    ifstream fin(full_file_name + manager_file_format);
        string buff;
        fin >> buff;
        if(buff == mark_file_manager_begin)
        {
            fin >> buff;
            if( buff == mark_file_commands)
                while( buff != mark_file_proceded_commands && !fin.eof())
                {
                    fin >> buff;    //cout << "BUFF\t" << buff << endl;
                    commands.push_back(buff);
                }
                commands.pop_back(); /// костыль
            if( buff == mark_file_proceded_commands)
                while( buff != mark_file_data_paths && !fin.eof())
                {
                    fin >> buff; //cout << "BUFF\t" << buff << endl;
                    proceded_commands.push_back(buff);
                }
                proceded_commands.pop_back(); /// костыль
            if( buff == mark_file_data_paths  && !fin.eof())
                while( buff != mark_file_manager_end)
                {
                    fin >> buff; //cout << "BUFF\t" << buff << endl;
                    data_paths.push_back(buff);
                }
                data_paths.pop_back(); /// костыль
        }
    fin.close();
}

void manager::Save_manager(string full_file_name)
{
    ofstream fout(full_file_name + manager_file_format);
        fout << mark_file_manager_begin << endl;
        fout << mark_file_commands << "\t";
            for (auto c : commands)
                fout << c << "\t";
                    fout << endl;
        fout << mark_file_proceded_commands << "\t";
            for (auto c : proceded_commands)
                fout << c << "\t";
                    fout << endl;
        fout << mark_file_data_paths << "\t";
            for (auto c : data_paths)
                fout << c << "\t";
                    fout << endl;
        fout << mark_file_manager_end;
    fout.close();
}

///EXTRA
string manager::get_info()
{
    string all_commands,all_proceded_commands,all_data_paths;
    for (auto c : commands)
            all_commands += (c + "\t");
    for (auto c : proceded_commands)
            all_proceded_commands += (c + "\t");
    for (auto c : data_paths)
            all_data_paths += (c + "\t");
    return "Manager object\t sucess=" + my_to_string(sucess) + "\t" + to_string(commands.size()) + "\t" + to_string(proceded_commands.size())
        + "\n\tcommands:\t" + all_commands + "\n\tproceded_commands:\t" + all_proceded_commands + "\n\tdata_paths:\t" + all_data_paths + "\n";
}
void manager::info()
{
    cout << "Manager object\t sucess=" << sucess << "\t" << commands.size() << "\t" << proceded_commands.size() << endl;
        cout << "\tcommands:\t";
        for (auto c : commands)
            cout << c << "\t";
        cout << "\n\tproceded_commands:\t";
        for (auto c : proceded_commands)
            cout << c << "\t";
        cout << "\n\tdata_paths:\t";
        for (auto c : data_paths)
            cout << c << "\t";
                cout << endl;
}

void manager::info_load_from_file()
{
    cout << "Manger load file info" << endl;
    cout <<  "\tbegins with \"" << mark_file_manager_begin << "\""<< endl;
    cout <<  "\tcommands load after \"" << mark_file_commands << "\"" << endl;
    cout <<  "\tproceded_commands load after \"" << mark_file_proceded_commands << "\"" << endl;
    cout <<  "\tends with \"" << mark_file_manager_end << "\"" << endl;
    cout << "WARNING data must end by symbol" << endl;
}

void manager::create_load_from_file_fishfile(const string full_file_name)
{
    ofstream fout(full_file_name + manager_file_format);
        fout << mark_file_manager_begin << endl;
        fout << mark_file_commands << "\t***enter_commands_instead***" << endl;
        fout << mark_file_proceded_commands << "\t***enter_proceded_commands_instead***" << endl;
        fout << mark_file_data_paths << "\t***enter_data_paths_instead***" << endl;
        fout << mark_file_manager_end << endl;
    fout.close();
}

///********************************************************************************************************

string manager_int_num(int i)
{//cout << "string int_num(int i)" << i << endl;
    string str_num = "";
    div_t d;    d = div(i, 10);
    if( d.quot > 0)
    {
        str_num = (int)d.rem+48;//cout << "str_num = " << str_num << endl;
        str_num = manager_int_num(d.quot) + str_num;
    }else{
        str_num = (int)d.rem+48;
         return str_num;
    }
    return str_num;
}

string my_to_string(double d)
{///cout << "FUNCTION my_to_string!!!!\n";
    if(int(d) == d)
        return manager_int_num(int(d));
    char arr_c [254]; sprintf(arr_c, "%f", d);
        string str_num = string(arr_c); return str_num;
}

string my_to_string(bool b)
{
    if(b)
        return "true";
    return "false";
}

double Load_marked_double_from_file(string Data_Way, string file_name, string mark, int index)
{//cout << "Load_marked_double_from_file()\n";
   // if(FileExists( (Data_Way + file_name) ){
        ifstream fin(Data_Way + file_name);
        string s_buff = "";
        double d_buff = 0;
        while(!fin.eof())
        {
            fin >> s_buff;//cout << mark << "\t" << s_buff << endl;
            if(s_buff == mark)
            {
                if(index != 0)
                    for(int i = 0; i < index; i++)
                        fin >> s_buff;
                fin >> d_buff; //cout << d_buff << endl;
                //break;
                fin.close();
                return d_buff;
            }
        }
        fin.close();
        //cout << d_buff << endl;
        return d_buff;
    /*}else{
        cout << "Load_marked_double_from_file no file: " << Data_Way + file_name << "\t0 returned\n";
        return 0;
    }*/
}

double Load_marked_double_from_file(string Data_Way, string file_name, string mark)///index = 1
    { return Load_marked_double_from_file(Data_Way, file_name, mark,0);}


///******************
///OSCILLATION

queue<manager_instruction> manager_read(string full_manager_file_name)
{
     ///очередь исполнения
    std::queue<manager_instruction> queue_instr;
    ///вектор потоков
    std::vector<std::thread> dynamic_thread;

    ///открытие инструкций к исполнению
    string s_buff;
    ifstream in_manager_instr(full_manager_file_name);
        while(s_buff != "START_INSTRUCTIONS" && !in_manager_instr.eof())
            in_manager_instr >> s_buff;///buff
        ///чтение имени файлов для обработки
        ///формирование очереди для процедуры обработки
        manager_instruction buff_instr;
        while(s_buff != "END_INSTRUCTIONS" && !in_manager_instr.eof())
        {
            in_manager_instr >> s_buff; //cout << s_buff << endl;
                buff_instr.summary_file = s_buff;
            in_manager_instr >> s_buff; //cout << s_buff << endl;
                buff_instr.DataWay = s_buff;
            in_manager_instr >> s_buff; //cout << s_buff << endl;
                buff_instr.base_file_name = s_buff;
            queue_instr.push(buff_instr);
        }
    in_manager_instr.close();//cout << endl << endl << endl;
    return queue_instr;
}

void manager_amplitude_analysis(queue<manager_instruction> queue_instr)
{cout << "manager_amplitude_analysis(queue<manager_instruction> queue_instr)" << endl;
    clock_t time_start = (double)clock()/CLOCKS_PER_SEC;
    try
    {
        if(queue_instr.front().base_file_name != "amplitude_analisys")
        {//cout << "throw\t" << queue_instr.front().base_file_name << endl; cerr << "manager wrong procedure !\n";
            throw "wrong_procedure\n";
        }
    }
    catch(char const* err)
    {
        cout << "wrong_procedure, process returned\n";
        return;
    }
    ///удаление проверочного элемента
    queue_instr.pop();
    ///вектор потоков
    std::vector<std::thread> dynamic_thread;/// исправить

    string base_file_name, Data_Way;
    ///managing loop
    while(queue_instr.size()>0 && queue_instr.front().base_file_name != "END_INSTRUCTIONS")
    {
        Data_Way = queue_instr.front().DataWay;
        base_file_name = queue_instr.front().base_file_name;
            //cout << DataWay << "\t" << base_file_name << endl;
            queue_instr.pop();
        dynamic_thread.push_back(std::thread(&General_Procedure_Oscillation,queue_instr.front().summary_file, Data_Way + base_file_name +"/", base_file_name, "_asb", ""));
        cout << "====================NEW TREAD====================\n\t\t\t" << dynamic_thread.back().get_id() << endl;
    }
        for(auto &th : dynamic_thread )
            th.join();
    cout << "manager_amplitude time used = " << (double)clock()/CLOCKS_PER_SEC - time_start << endl;
};

std::mutex mtx;
void General_Procedure_Oscillation(string DataWay_summary, string Data_Way, string base_file_name, string name_add, string mode_middle_angle_correction)
{
    std::lock_guard<std::mutex> lg(mtx);
    oscillation D;
    matrix b_angle, old_a;
    int index_balance_column = 4;
    ///загрузка подготовленных данных колебаний (без потока)
    D.Load_Oscillation_Prop_File(Data_Way, base_file_name + "_r");
    ///загрузка идентификатора данных
    D.set_data_source(base_file_name);
    ///загрузка готовых данных потока (по идентификатору данных)
    //D.Load_flow_summary(DataWay_summary, "overall_flow_summary");///пересмотреть, нужно оптимизировать
    ///загрузка готовых данных модели (по модели из данных потока)
    D.Load_model_info(DataWay_summary, "model_info");

        ///процедура анализа амплитуды
        procedure_ampl_analisys(Data_Way, base_file_name + "_r", "", D, 25.0, 0.0);

    if( mode_middle_angle_correction != "")
    {
      ///загрузка балансировочного угла
        b_angle.Load_matrix(Data_Way + base_file_name + "_r" + "_middle_angle" + ".txt");
        D.set_balance_angle(b_angle);
        ///смещение угла
        D.angle_shift_balance(index_balance_column);
        ///процедура анализа амплитуды исправленных углов
        procedure_ampl_analisys(Data_Way, base_file_name + "_r", name_add, D, 25.0, 0.0);
        ///сохранение результатов
        D.Save_Oscillation(Data_Way, base_file_name + name_add + "_r");
    }

        ///кат данных по инструкции
        //Redd_Cut_file(Data_Way, "Proc_cut");
        ///создание инструкции для Gnuplot
        //Make_Gnuplot_File_Pic_Plotter(base_file_name, "", Data_Way + "Pic_plot");
}


///Fluent report
void procedure_get_my_res(string Data_Way,string Data_Way_Res, string base_file_name, double angle)
{
    ifstream fin(Data_Way + base_file_name + ".txt");
    cout << Data_Way + base_file_name + ".txt\n";
    string buff;
    double cx_all, cx_aft1, cx_aft2, cx_front1, cx_front2, cy_all, cy_aft1, cy_aft2, cy_front1, cy_front2, mz_all, mz_aft1, mz_aft2, mz_front1, mz_front2; ///из-за скобочки тк вид следующий: 0.034324234)

    ///вид следующий:
    ///Force vector: (1 0 0)
    ///body_aft_2   6_столбцов
    ///body_front_2 6_столбцов
    ///body_aft_1   6_столбцов
    ///body_front_1 6_столбцов
    fin >> buff;cout << buff << endl;
    //cout << "Force vector: (1 0 0)\n";
    while(buff != "body_aft_2")
        {fin >> buff;/*cout << buff << endl;*/}cout << "buff" << endl;
            fin >> buff >> buff >> buff >> buff >> buff >> cx_aft2;
    while(buff != "body_front_2")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cx_front2;
    while(buff != "body_aft_1")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cx_aft1;
    while(buff != "body_front_1")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cx_front1;
    while(buff != "net")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cx_all;

    cout << "Force vector: (0 1 0)\n";
    ///Force vector: (0 1 0)
    while(buff != "body_aft_2")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cy_aft2;
    while(buff != "body_front_2")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cy_front2;
    while(buff != "body_aft_1")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cy_aft1;
    while(buff != "body_front_1")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cy_front1;
    while(buff != "net")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> cy_all;

    ///Moment Center: (0.091499999 0 0)
    while(buff != "body_aft_2")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >> buff >> buff >> buff;
        buff.erase(buff.find(')'), 1);
        mz_aft2 = strtod(buff.c_str(), NULL);

    while(buff != "body_front_2")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >> buff >> buff >> buff;
        buff.erase(buff.find(')'), 1);
        mz_front2 = strtod(buff.c_str(), NULL);

    while(buff != "body_aft_1")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >> buff >> buff >> buff;
        buff.erase(buff.find(')'), 1);
        mz_aft1 = strtod(buff.c_str(), NULL);

    while(buff != "body_front_1")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >> buff >> buff >> buff;
        buff.erase(buff.find(')'), 1);
        mz_front1 = strtod(buff.c_str(), NULL);

    while(buff != "net")
        fin >> buff;
            fin >> buff >> buff >> buff >> buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >>buff >> buff >> buff >> buff >> buff >> buff;
        buff.erase(buff.find(')'), 1);
        mz_all = strtod(buff.c_str(), NULL);

    fin.close();

    ///запись
    ofstream fout(Data_Way_Res + "MyRes_" + ".txt", ios::app);
        fout << angle << "\t" << 0 << "\t\t" <<
            cx_all << "\t" << cy_all << "\t" << mz_all << "\t\t" <<
                cx_aft1 << "\t" << cy_aft1 << "\t" << mz_aft1 << "\t\t" <<
                    cx_aft2 << "\t" << cy_aft2 << "\t" << mz_aft2 << "\t\t" <<
                        cx_front1 << "\t" << cy_front1 << "\t" << mz_front1 << "\t\t" <<
                            cx_front2 << "\t" << cy_front2 << "\t" << mz_front2 << "\t\t" <<
                                (double)mz_aft1 + (double)mz_aft2 << endl;
    fout.close();
};

void symmetry_coeff_report_flu(string Data_Way, string file_name)
{
    matrix D;
        D.Load_matrix(Data_Way + file_name + ".txt");
    for(int i = 0; i< D.get_height(); i++)
    {
        ///angle
        D.set_data(i, 0, D.get_data(i, 0)*(-1));
        ///aft1 cy mz
        D.set_data(i, 3, D.get_data(i, 3)*(-1));
        D.set_data(i, 4, D.get_data(i, 4)*(-1));
        ///aft2 cy mz
        D.set_data(i, 6, D.get_data(i, 6)*(-1));
        D.set_data(i, 7, D.get_data(i, 7)*(-1));
        ///front1 cy mz
        D.set_data(i, 9, D.get_data(i, 9)*(-1));
        D.set_data(i, 10, D.get_data(i, 10)*(-1));
        ///front2 cy mz
        D.set_data(i, 12, D.get_data(i, 12)*(-1));
        D.set_data(i, 13, D.get_data(i, 13)*(-1));
        ///aft1 cy mz
        D.set_data(i, 15, D.get_data(i, 15)*(-1));
        D.set_data(i, 17, D.get_data(i, 17)*(-1));
        ///aft_sum mz
        D.set_data(i, 18, D.get_data(i, 18)*(-1));
    }
    D.Save_matrix(Data_Way + "sym_"+ file_name + ".txt");
}

void Flow_pulsation_check(string Data_Way, string file_name)
{
    matrix A;///iter / data
        A.Load_matrix(Data_Way + file_name);
    double avg = 0;
        for(int i = 0; i < A.get_height(); i++)
            avg += A.get_data(i,1);
                avg = avg/A.get_height();
                    cout << avg << endl;

    for(int i = 0; i < A.get_height(); i++)
    {
        A.set_data(i,0, i);
        A.set_data(i,1, A.get_data(i,1)-avg);
    }
    avg = 0;
    A.info();
        A.Save_matrix(Data_Way + file_name + "_f.txt");

    matrix W;
    //Get_Spectr(A, W, A.get_height(), 4);

    W = Get_Spectr_column_gap(A.get_column(0), A.get_column(1), 0, A.get_height(), 4);
    for(int i = 0; i < W.get_height(); i++)
       W.set_data(i,1, log(W.get_data(i,1)));
       //W.set_data(i,1, W.get_data(i,1));
    W = W.get_matrix_part(0,W.get_height()/2, 0,2);
    W.Save_matrix(Data_Way + file_name + "_spec.txt");
}

///INERTIA

///Поправка случаев когда датчик зашкаливает и переходит на следующий оборот
void Fix_LIR_Data(string Data_Way, string base_file_name, string result_file, double border)
{
    matrix A1, A2;
        A1.Load_matrix(Data_Way + base_file_name + ".txt");
        A2 = A1;
    double data_check = 0;
    bool flag = false; ///нужно ли выполнять смещение участка?

    for(int i = 1; i < A1.get_height(); i++)
    {cout << i << endl;
        data_check = A1.get_data(i-1,1) - A1.get_data(i,1);
        if(data_check < -1.0*border )
            flag = true;
        if(data_check > border )
            flag = false;
        if(flag==true )
            A2.set_data(i,1, A1.get_data(i,1) - 360.0);
    }
    A2.info();
    A2 = A2.m_column_shift(1, -A2.get_data(A2.get_height()-1, 1));
    A2.Save_matrix(Data_Way + result_file + ".txt");
}

void procedure_inertia(string Data_Way, string base_file_name, int  N_p,  int overlap, int fft_zero_index, double m1, double r1, double g)
{
    oscillation Inert1;
        Inert1.Load_Oscillation_Protocol(Data_Way, base_file_name);
            Inert1.Save_Oscillation(Data_Way, base_file_name + "_osc.txt");
        matrix in1; in1 = (Inert1.get_time()).merge_width(Inert1.get_angle(), 1);in1.info();
        matrix S, S_m, data_calc(1,3);  ///data_calc() amplitude, I1, I2
        ///амплитуду надо получать из огибающей.

        double T0 = 0, amplitude = 0;
        int index_h = 0, step = (int)Inert1.get_angle().get_height()/N_p;
        for(int i = 0; i < N_p; i++)
        {
            S = Get_Spectr_column_gap(Inert1.get_time(), Inert1.get_angle(), index_h, index_h + step, fft_zero_index);
                S.Save_matrix(Data_Way + base_file_name + "_fft" + to_string(i) + ".txt");
            ///время/частота
            S_m = Get_Spectr_max_w(S);
                S_m.set_data(0,0, Inert1.get_time().get_data(index_h,0));

                    T0 = 1.0/S_m.get_data(0,1);
                    amplitude = 53.3-1.3*S_m.get_data(0,0);///амплитуду надо получать из огибающей. ///(Inert1.get_angle()).get_data( (index_h + (int)step/2), 0);

                    data_calc.set_data(0,0, amplitude);///amplitude
                        data_calc.set_data(0,1, r1*m1*g*pow(T0,2)/4.0/pow(PI,2));///I1
                            data_calc.set_data(0,2, data_calc.get_data(0,1)/(pow((1+pow(sin(amplitude/2.0*0.0174533), 2)/4.0),2)) );///I2

                    ///присоединение справа столбцов счета
                    S_m = S_m.merge_width(data_calc, S_m.get_width());
                        S_m.Save_matrix_end(Data_Way + base_file_name +"_res.txt");
            index_h += step;
        }
}

void procedure_inertia(string Data_Way, string base_file_name, int  N_p, int overlap, int fft_zero_index, double sigma, int K_pow, double m1, double r1, double g)
{
    oscillation Inert1;
        Inert1.Load_Oscillation_Protocol(Data_Way, base_file_name);
            Inert1.Save_Oscillation(Data_Way, base_file_name + "_osc.txt");//matrix in1; in1 = (Inert1.get_time()).merge_width(Inert1.get_angle(), 1);in1.info();
        matrix S, S_m, Res, data_calc(1,3);  ///data_calc() amplitude, I1, I2
        ///амплитуду надо получать из огибающей.
        matrix ampl_env_top, ampl_env_bot, dampl;
            dampl = derevative_3_dot(Inert1.get_angle(), Inert1.get_time());
                ampl_env_top = Get_Envelope_top( (Inert1.get_time()).merge_width(Inert1.get_angle(), (Inert1.get_time()).get_height()), dampl);
                ampl_env_bot = Get_Envelope_bot( (Inert1.get_time()).merge_width(Inert1.get_angle(), (Inert1.get_time()).get_height()), dampl);
                    ampl_env_top.Save_matrix(Data_Way + base_file_name + "_ampl_envelop_top.txt");
                    ampl_env_bot.Save_matrix(Data_Way + base_file_name + "_ampl_envelop_bot.txt");
                    ///строить линейную аппроксимацию для описания амплитуды

        double T0 = 0, amplitude = 0, K = 0, Iavg = 0;
        int index_h = 0, window = (int)Inert1.get_angle().get_height()/N_p, I_counter = 0;
        //int window_mode = 1, window_counter = 0;

        cout << "window = " << window << "\tstep =" << window-overlap << endl << endl;

        string get_freq_way;
            get_freq_way = "front";
            //get_freq_way = "FFT";
        double min_ampl = 23.0, incline_curr = 0;
        //int signal_front_counter = 0;
        int i_curr = 0; /// текущий отсчет огибающей
        while( index_h + window-overlap < Inert1.get_angle().get_height())
        { cout << "Procedure step\t time = " << Inert1.get_time().get_data(index_h,0) << endl;
            T0 = 0;
            ///период через спектры
            if(get_freq_way == "FFT")
            {
                S = Get_Spectr_column_gap(Inert1.get_time(), Inert1.get_angle(), index_h, index_h + window, fft_zero_index, sigma);///S.Save_matrix(Data_Way + base_file_name + "_fft" + my_to_string(i) + ".txt");///сохранение спектров в файл
                S_m = Get_Spectr_max_w(S);
                S_m.set_data(0,0, Inert1.get_time().get_data(index_h,0));
                    T0 = 1.0/S_m.get_data(0,1);
            }

            ///период по фронтам
            if(get_freq_way == "front")
            {//cout << "index_h = " << index_h << endl;
                S = (Inert1.get_time()).merge_width(Inert1.get_angle(), 1);
                    S = S.get_matrix_part(index_h, index_h + window, 0,1); S.info();

                S_m = Front_Get_Main_freq(S, ampl_env_top, ampl_env_bot);
                    S_m.set_data(0,0, Inert1.get_time().get_data(index_h,0));
                    T0 = 1.0/S_m.get_data(0,1);
            }
            //cout << "Period founded from\t = " << S_m << endl;

            ///амплитуду надо получать из огибающей. ///(Inert1.get_angle()).get_data( (index_h + (int)window/2), 0);
            for(int i = 0; i < ampl_env_top.get_height(); i++)
                if(ampl_env_top.get_data(i,0) >= S_m.get_data(0,0))
                {//cout <<"i = " << i << "\t" << ampl_env_top.get_data(i,0) << "\t >=?\t" << S_m.get_data(0,0) << endl;
                    i_curr = i;
                    i = ampl_env_top.get_height();
                }
            Find_Incline_Coefficient(ampl_env_top.get_matrix_part(i_curr, i_curr+2, 0, 2), &incline_curr, 2);
                //cout << "Incline found from matrix = \n" << ampl_env_top.get_matrix_part(i_curr, i_curr+2, 0, 2) << "\tincline = " << incline_curr << endl;
            amplitude = ampl_env_top.get_data(i_curr, 1) ;//+ incline_curr*S_m.get_data(0,0);
               // cout << "ampl0 = " << ampl_env_top.get_data(i_curr, 1) << "\tampl1 = " << ampl_env_top.get_data(i_curr+1, 1) << "\ttime = " << S_m.get_data(0,0) << "\tamplitude_current= " <<  amplitude << endl;
           // cout << "*****************************\n\n\n";
            ///amplitude = 52.5-1.4*S_m.get_data(0,0);
                    ///полный нормальный эллиптический интеграл Лежандра 1-го рода до порядка
                    K = Legander_K(amplitude/2.0*0.0174533, K_pow);

                    data_calc.set_data(0,0, amplitude);///amplitude
                        data_calc.set_data(0,1, r1*m1*g*pow(T0,2)/4.0/pow(PI,2));///I1
                            data_calc.set_data(0,2, data_calc.get_data(0,1)/pow(K,2) );///I2
                                if(amplitude > min_ampl)
                                {
                                    Iavg+=data_calc.get_data(0,2);
                                    I_counter++;
                                }
                ///присоединение справа столбцов счета
                S_m = S_m.merge_width(data_calc, S_m.get_width());///S_m.Save_matrix_end(Data_Way + base_file_name +"_res.txt");
                    Res = Res.merge_height(S_m, Res.get_height());
        index_h += (window - overlap);
        }
        ///огибающая полученного момента инерции
        matrix Env_I;
        matrix dRes; dRes = derevative_3_dot(Res, 4, 0);
            Env_I = Get_Envelope_top(Res, 4, 0, dRes);
                Env_I.Save_matrix(Data_Way + base_file_name + "_I_envelop_top.txt");
            Env_I = Get_Envelope_bot(Res, 4, 0, dRes);
                Env_I.Save_matrix(Data_Way + base_file_name + "_I_envelop_bot.txt");

        Res.Save_matrix(Data_Way + base_file_name +"_res.txt");

        matrix Overall_Res(1,5);
            Overall_Res.set_data(0,0, window);
            Overall_Res.set_data(0,1, overlap);
            Overall_Res.set_data(0,2, K_pow);
            Overall_Res.set_data(0,3, min_ampl);
            Overall_Res.set_data(0,4, Iavg/I_counter);
                Overall_Res.Save_matrix_end(Data_Way + base_file_name + "_overall_res.txt");
        cout << "summ = " << Iavg << "\tI_counter = " << I_counter << "\tIavg = " << Iavg/I_counter << endl;
}

matrix procedure_get_pendulum_frequency(string Data_Way, string base_file_name, oscillation D, string method_find_freq, int window_size, int step_size, int zero_index, int K_power, double avg_angle)///time // freq ///f = 1/T
{cout << "angle_changing\t" << avg_angle << endl;
    ///вычет угла
    matrix a; a.info(); a = D.get_angle();
        for(int i = 0; i < a.get_height(); i++)
            a.set_data(i,0, a.get_data(i,0)-avg_angle);
        //a.m_column_shift(0, -1.0*avg_angle);
            D.set_angle(a);
    return procedure_get_pendulum_frequency(Data_Way,  base_file_name, D,  method_find_freq, window_size, step_size, zero_index, K_power);///time // freq ///f = 1/
}

/*matrix procdure_get_pendulum_frequency(string Data_Way, string base_file_name, oscillation D, string method_find_freq, int window_size, int step_size, int zero_index, int K_power)///time // freq ///f = 1/T
{cout << "get_avg_angle()\n";
  ///вычетание средней сигнала
    double avg_angle = 0;
    for(int i = 0; i < (D.get_angle()).get_height(); i++)
    {
        avg_angle += D.get_angle().get_data(i,0);
    }
        avg_angle = avg_angle/(D.get_angle()).get_height();
                    cout << "FFT: avg_angle = " << avg_angle << endl;
                        ofstream fout(Data_Way + base_file_name + "_w_avg_angle.txt"); fout << avg_angle; fout.close();
    return procdure_get_pendulum_frequency(Data_Way,  base_file_name, D,  method_find_freq, 0.0, window_size, step_size, zero_index, K_power);///time // freq ///f = 1/T
}*/

matrix procedure_get_pendulum_frequency(string Data_Way, string base_file_name, oscillation D, string method_find_freq, int window_size, int step_size, int zero_index, int K_power)///time // freq ///f = 1/T
{cout << "procedure_get_pendulum_frequency()\n";

    string pendulum_name = "_pendulum";
    matrix f_m, dat; ///рузльтатирующая матрица время/частота
    ///Нахождение частоты реализации
    ///Реализация делится на окна указанной ширины window_size, кроме того окно движется с шагом step_size
    //cout << "Save_Oscillation\n";
    D.get_angle().Save_matrix(Data_Way + base_file_name + "_f1.txt");

    if(/*!FileExists(Data_Way + base_file_name + "_w0.txt")*/true)
    {
        if (method_find_freq == "FFT")
        {
            dat = D.get_time();
            dat = dat.merge_width(D.get_angle(), 1);
            f_m = Get_Freq_Window_FFT(dat, step_size, window_size, zero_index);
                /**
                matrix spectr, f;
            //int window_counter = 0;
            for(int index_h = 0; index_h < D.get_angle().get_height()-window_size;)
            {cout << "FFT index_h = " << index_h << "\t" << (double)index_h/D.get_angle().get_height()*100 << "%\t";
                spectr = Get_Spectr_column_gap(D.get_time(), D.get_angle(), index_h, index_h + window_size, zero_index);
                    f = Get_Spectr_max_w(spectr);///плотность и наиболее плотная частота между index_h и index_h + window_size
                    f.set_data(0,0, D.get_time().get_data(index_h + (int)window_size/2,0));///время в середине участка и частота главной гармоники ///убрано!!!!!!!
                f_m = f_m.merge_height(f, f_m.get_height());///собираем полученные частоты
                index_h +=step_size;
                cout  << "\tfinished\n";
            }
                */
        }//f_m.info();
    f_m.Save_matrix(Data_Way + base_file_name + "_w0.txt");
    }else{
        f_m.Load_matrix(Data_Way + base_file_name + "_w0.txt");
    }

    cout << "FFT finished\n";
    /**if(method == "front")
    {
        return 0;
    }*/

    matrix dangle; dangle = derevative_3_dot(D.get_angle(), D.get_time());//dangle.info();
        matrix angle_envelop_top;
            angle_envelop_top = Get_Envelope_top(D.get_angle(), D.get_time(), dangle.get_column(1));
                    //cout << "envelop_top\n";angle_envelop_top.info();
                angle_envelop_top.Save_matrix(Data_Way + base_file_name + "_w_envelop_top.txt");
        /**matrix angle_envelop_bot;
            angle_envelop_bot = Get_Envelope_top((D.get_time()).merge_width(D.get_angle(), (D.get_angle()).get_height()), dangle);*/

    ///получение амплитуды от вермени(!) через огибающую и линейную аппроксимацию между точек
    matrix amplitude_line(angle_envelop_top.get_height()-1, 2);///t[i] t[i+1] a b ///a*x+b
    double a = 0, b = 0;
    matrix angle_envelop_part;
    double arr_a [amplitude_line.get_height()];
    double arr_b [amplitude_line.get_height()];
    for(int i = 0; i < amplitude_line.get_height(); i++)///строим аппроксимации между точек огибающей
    {//cout << "amplitude_line i = " << i << "\t";
        amplitude_line.set_data(i, 0, angle_envelop_top.get_data(i,0));///time[i]
            amplitude_line.set_data(i, 1, angle_envelop_top.get_data(i+1,0));///time[i]
        angle_envelop_part = angle_envelop_top.get_matrix_part(i, i+2, 0, 2);
            Find_Incline_Shift_Coefficient(angle_envelop_part, &a, &b, 2);
        //amplitude_line.set_data(i, 2, a);///a[i]
          //  amplitude_line.set_data(i, 3, b);///b[i]
        arr_a[i] = a;
        arr_b[i] = angle_envelop_top.get_data(i,1);//b;
    }//cout << "line approx finished\n";
    a = 0; b = 0;

    matrix n_h;/// коэффициенты a и b
    Filter_Outlier(arr_a, 0, amplitude_line.get_height(), 1.4, 2.5);
        n_h.array_to_matrix_column(arr_a,amplitude_line.get_height());
            amplitude_line = amplitude_line.merge_width(n_h, amplitude_line.get_width());
    Filter_Outlier(arr_b, 0, amplitude_line.get_height(), 1.4, 2.5);
        n_h.array_to_matrix_column(arr_b,amplitude_line.get_height());
            amplitude_line = amplitude_line.merge_width(n_h, amplitude_line.get_width());
    //cout << "amplitude line\n";
    //amplitude_line.info();
    amplitude_line.Save_matrix(Data_Way + base_file_name + "_w_amplitude_line.txt");
    //cout << "Filtering Outlier finished\n";

    matrix amplitude(f_m.get_height(),1);///доп инфа к f_m
    double t0 = 0;
    a = amplitude_line.get_data(0,2);
    b = amplitude_line.get_data(0,3);
        t0 = amplitude_line.get_data(0,1)-amplitude_line.get_data(0,0);
    amplitude.set_data(0,0, a*(f_m.get_data(0,0)-t0)+b);
    //cout << "amplitude counted \n";


    //cout << "i = " << 0 <<  "\ta*t+b " << a << "*" << (f_m.get_data(0,0)-t0) << " + " << b << " = " << a*(f_m.get_data(0,0)-t0)+b << endl;
    for(int i = 1; i < amplitude.get_height(); i++)
    {
        for(int k = 0; k < amplitude_line.get_height()-1; k++)///ищем нужный промежуток
        {
            if(amplitude_line.get_data(k,0) <= f_m.get_data(i,0) && amplitude_line.get_data(k, 1) >= f_m.get_data(i,0) )
            {
                a = amplitude_line.get_data(k,2);
                b = amplitude_line.get_data(k,3);
                t0 = amplitude_line.get_data(k,0);
                k = amplitude_line.get_height();
                //cout << a << "\t" << b << endl;
            }
        }
        amplitude.set_data(i,0, a*(f_m.get_data(i,0)-t0)+b);
            //cout << "i = " << i <<"\tt0=" << t0 <<  "\tt = " << f_m.get_data(i,0) << "\ta*t+b " << a << "*" << f_m.get_data(i,0)-t0 << " + " << b << " = " << a*(f_m.get_data(i,0)-t0)+b << endl;
        //cout << "amplitude.get_data(i,0) = " << amplitude.get_data(i,0) << endl;
    }
    ///поправка к частоте для физического маятника (полином Л)
    matrix fixed_freq(f_m.get_height(), 1);
    for(int i = 0; i < fixed_freq.get_height(); i++)
        fixed_freq.set_data(i, 0, f_m.get_data(i, 1)*Legander_K(amplitude.get_data(i,0)/2.0*0.0174533, K_power));
    f_m = f_m.merge_width(fixed_freq, f_m.get_width());
    f_m = f_m.merge_width(amplitude, f_m.get_width());
        f_m.Save_matrix(Data_Way + base_file_name + pendulum_name + ".txt");
    return f_m;
}


void procedure_moment_of_inertia_stand(string Data_Way, string base_file_name, string method_find_freq, int window_size, int step_size,int zero_index, int K_power)
{
    oscillation D;
        D.Load_Oscillation_Protocol(Data_Way,  base_file_name);
    matrix pendulum_data;///time/freq_mesured/freq_fixed/amplitude
        string pendulum_name = "_pendulum";
        pendulum_data.Load_matrix(Data_Way + base_file_name + pendulum_name + ".txt");
            ///при отсутсвии идет счет, при наличии загрузка.
            /*    if(!FileExists(Data_Way + base_file_name + ".txt"))
            {cout << "Count Pendulum\n";
                pendulum_data = procdure_get_pendulum_frequency(Data_Way, base_file_name, D, method_find_freq, window_size, step_size, zero_index, K_power);
            }else{cout << "Load Pendulum\n";
                pendulum_data.Load_matrix(Data_Way + base_file_name + pendulum_name + ".txt");
            }*/
    ///averaging -> freq
        ///creteria
        /// ampl_min <= amplitude <= ampl_max
    double ampl_max = 60, ampl_min = 8;
    ///**********************

    int h = 0, i = 0, i0 = 0;
    while(pendulum_data.get_data(i,3) >= ampl_max && i < pendulum_data.get_height()-1)
        i++;
    i0 = i;
    while(pendulum_data.get_data(i,3) >= ampl_min && i < pendulum_data.get_height()-1)
    {//cout << i << endl;
        i++;
        h++;
    }
    cout << "i0 = " << i0 << "\th = " << h << endl;
    matrix work_freq(h,2);
    matrix work_ampl(h,1);
    ///double* ampl_arr = new double [h];
    for(int k = 0; k < h; k++)
    {
        work_freq.set_data(k, 0, pendulum_data.get_data(k+i0, 0));
        ///ampl_arr[k] = pendulum_data.get_data(k+i0,2);
        work_freq.set_data(k, 1, pendulum_data.get_data(k+i0, 2));
        work_ampl.set_data(k, 0, pendulum_data.get_data(k+i0, 3));
    }
    ///Filter_Outlier1(ampl_arr, h, 0.05, 0.05);
    ///matrix ampl_column;
       /// ampl_column.array_to_matrix_column(ampl_arr, h); ampl_column.info();
    ///work_freq = work_freq.merge_width(ampl_column, work_freq.get_width());
    work_freq = smooth_matrix_column(work_freq, 10, 1);
        work_freq = smooth_matrix_column(work_freq, 10, 1);
            work_freq = smooth_matrix_column(work_freq, 10, 1);
    work_freq = work_freq.merge_width(work_ampl, work_freq.get_width());
    work_freq.Save_matrix(Data_Way + base_file_name + "_filt.txt");

    ///получение частоты f0
   /* while()


    ///properties
    ///m g r
    double m = 1.375, g = 9.83202, r = 0.00423139, f0;

    ///result
    double I;
    I = m*g*r/4/PI/PI/f0/f0;*/
}

matrix procedure_osc_analisys(string Data_Way, string base_file_name, oscillation D, int window_size, int step_size, int zero_index)///08.01.20
{
    ofstream fout_input(Data_Way + base_file_name + "_input.txt");
        fout_input << Data_Way + base_file_name <<  "\nwindow_size = " << window_size << "\nstep_size = " << step_size << "\nzero_index = " << zero_index;
    fout_input.close();
    matrix f_m, dat; ///рузльтатирующая матрица время/частота
    ///Нахождение частоты реализации
    ///Реализация делется на окна указанной ширины window_size, кроме того окно движется с шагом step_size
    //cout << "Save_Oscillation\n";
    D.get_angle().Save_matrix(Data_Way + base_file_name + "_f1.txt");

    if(/*!FileExists(Data_Way + base_file_name + "_w0.txt")*/true)
    {
        //if (method_find_freq == "FFT")
        //{
        ///фурье от сигнала
            dat = D.get_time();
            dat = dat.merge_width(D.get_angle(), 1);
                f_m = Get_Freq_Window_FFT(dat, step_size, window_size, zero_index);
        //}//f_m.info();
                    f_m.Save_matrix(Data_Way + base_file_name + "_w0.txt");

        ///фурье от производной сигнала
            dat = D.get_time();
            dat = dat.merge_width(D.get_dangle(), 1);
                f_m = Get_Freq_Window_FFT(dat, step_size, window_size, zero_index);
                    f_m.Save_matrix(Data_Way + base_file_name + "_d1w0.txt");

    }else{
        cout << "File Found! Loading Freq-file...\n";
        f_m.Load_matrix(Data_Way + base_file_name + "_w0.txt");
    }
    matrix m; return m;
}

matrix procedure_Redd(matrix envelop, oscillation D)
{
    matrix data_Redd(envelop.get_height()-2,3);
    bool key1 = false, key2 = false;
    ///проверка определенности
    if(D.get_flow_v(0)!=0)
    {
        key1 = true;
        cout << "Flow is defined\n";
    }else{
        key1 = false;
        cout << "Flow is indefined\n";
    }
    if(D.get_model_S()!=0)
    {
        key2 = true;
        cout << "Model parameters is defined\n";
    }else{
        key2 = false;
        cout << "Model parameters is indefined\n";
    }
    if(key1 && key2)
        cout << "m_dyn coeff = " << D.get_flow_v(0)/D.get_flow_q(0)*D.get_model_I()/D.get_model_S()/D.get_model_l()/D.get_model_l() * 4.0 << endl;

        //cout<< key1 << key2 << endl;
        //cout << D.get_flow_v(0) << "\t" << D.get_flow_q(0) << "\t" << D.get_model_I() << "\t" << D.get_model_S() << "\t" << D.get_model_l() << endl;
        cout << "m_dyn coeff = " << D.get_flow_v(0)/D.get_flow_q(0)*D.get_model_I()/D.get_model_S()/D.get_model_l()/D.get_model_l() * 4.0 << endl;

    for(int i = 0; i < data_Redd.get_height(); i++)
    {
        data_Redd.set_data(i,0, envelop.get_data(i,0));///time
        data_Redd.set_data(i,1, envelop.get_data(i,1));///angle
        if (envelop.get_data(i+1,1)/envelop.get_data(i,1) > 0.0)
        {
            if( !key1 && !key2)///flow indef && model indef
                data_Redd.set_data(i,2, 2.0*log(envelop.get_data(i+1,1)/envelop.get_data(i,1))/(envelop.get_data(i+1,0)-envelop.get_data(i,0)));
            if( key1 && !key2)///flow def && model indef *v/q
                data_Redd.set_data(i,2, D.get_flow_v(i)/D.get_flow_q(i) * 2.0*log(envelop.get_data(i+1,1)/envelop.get_data(i,1))/(envelop.get_data(i+1,0)-envelop.get_data(i,0)));
            if( !key1 && key2)///flow indef && model def *I/S/l/l
                data_Redd.set_data(i,2, D.get_model_I()/D.get_model_S()/D.get_model_l()/D.get_model_l() * 2.0*log(envelop.get_data(i+1,1)/envelop.get_data(i,1))/(envelop.get_data(i+1,0)-envelop.get_data(i,0)));
            if( key1 && key2)///flow def && model def
                data_Redd.set_data(i,2, D.get_flow_v(i)/D.get_flow_q(i)*D.get_model_I()/D.get_model_S()/D.get_model_l()/D.get_model_l() * 2.0*log(envelop.get_data(i+1,1)/envelop.get_data(i,1))/(envelop.get_data(i+1,0)-envelop.get_data(i,0)));
        }
        else
            data_Redd.set_data(i,2,0);
    }
    return data_Redd;
}

void procedure_ampl_analisys(string Data_Way, string base_file_name, string name_add, oscillation &D, double freq, double angle_shift)
{cout << "procedure_ampl_analisys(string Data_Way, string base_file_name, oscillation D, double freq)\n";
    base_file_name += name_add;

    ///ANGLE SHIFT
    if(angle_shift != 0.0)
    {cout << "Angle shift\n";
        for(int i = 0; i < D.get_angle().get_height()-1; i++)
            D.set_angle_data(i, 0, D.get_angle().get_data(i,0) + angle_shift);
    }
    ///ENVELOP
    cout<< "Envelop...\n";
    matrix env_t, env_t3, env_t10, env_b, env_b3, env_b10;
        env_t = Get_Envelop_frq_top(D.get_time(), D.get_angle(), freq, 0.5);
            env_t.Save_matrix(Data_Way + base_file_name + "_env_top.txt");
                env_t3 = smooth_matrix_column(env_t, 3, 1);
                    env_t3.Save_matrix(Data_Way + base_file_name + "_env_top_sm3.txt");
                env_t10 = smooth_matrix_column(env_t, 10, 1);
                    env_t10.Save_matrix(Data_Way + base_file_name + "_env_top_sm10.txt");

        env_b = Get_Envelop_frq_bot(D.get_time(), D.get_angle(), freq, 0.5);
            env_b.Save_matrix(Data_Way + base_file_name + "_env_bot.txt");
                env_b3 = smooth_matrix_column(env_b, 3, 1);
                    env_b3.Save_matrix(Data_Way + base_file_name + "_env_bot_sm3.txt");
                env_b10 = smooth_matrix_column(env_b, 10, 1);
                    env_b10.Save_matrix(Data_Way + base_file_name + "_env_bot_sm10.txt");
    ///единый записаный последоватльно по времени
    matrix env_tb, env_t3b3, env_t10b10;
        env_tb = Get_Envelope_topbot(0, env_t, env_b);
        env_t3b3 = Get_Envelope_topbot(0, env_t3, env_b3);
        env_t10b10 = Get_Envelope_topbot(0, env_t10, env_b10);
            env_tb.Save_matrix(Data_Way + base_file_name + "_env_tb.txt");
            env_t3b3.Save_matrix(Data_Way + base_file_name + "_env_tb_sm3.txt");
            env_t10b10.Save_matrix(Data_Way + base_file_name + "_env_tb_sm10.txt");

    ///анализ среднего угла 1-time, 2-avg_mid осреднено по всем углам, 3-env_mid осреднено по огибающим
    cout << "MIDDLE ANGLE PROCEDURE\n";
    matrix middle_angle1, middle_angle_sm1, middle_angle_sm2, middle_angle_sm3;
    int smooth_window_length = 1;
        D.middle_angle_timeavg_halfT(env_t, env_b,middle_angle1);
    cout << "MIDDLE ANGLE FINISHED\n";
    smooth_window_length = 4;
        middle_angle_sm1 = smooth_matrix_column(middle_angle1, smooth_window_length, 1);
            middle_angle_sm1.erase_column(0);
    smooth_window_length = 16;
        middle_angle_sm2 = smooth_matrix_column(middle_angle1, smooth_window_length, 1);
            middle_angle_sm2.erase_column(0);
    smooth_window_length = 64;
        middle_angle_sm3 = smooth_matrix_column(middle_angle1, smooth_window_length, 1);
            middle_angle_sm3.erase_column(0);
    middle_angle1 = middle_angle1.merge_width(middle_angle_sm1);
    middle_angle1 = middle_angle1.merge_width(middle_angle_sm2);
    middle_angle1 = middle_angle1.merge_width(middle_angle_sm3);

            //middle_angle_sm1.Save_matrix(Data_Way + base_file_name + "_test1.txt");
            //middle_angle_sm2.Save_matrix(Data_Way + base_file_name + "_test2.txt");
            //middle_angle_sm3.Save_matrix(Data_Way + base_file_name + "_test3.txt");
            middle_angle1.Save_matrix(Data_Way + base_file_name + "_middle_angle.txt");

    ///TEST
        ///top sm0
        procedure_Redd(env_t,D).Save_matrix(Data_Way + base_file_name + "_env_top_redd.txt");//cout << "Redd\tt0\n";
        ///top sm3
        procedure_Redd(env_t3,D).Save_matrix(Data_Way + base_file_name + "_env_top_sm3_redd.txt");//cout << "Redd\tt3\n";
        ///top sm10
        procedure_Redd(env_t10,D).Save_matrix(Data_Way + base_file_name + "_env_top_sm10_redd.txt");//cout << "Redd\tt10\n";
        ///bot sm0
        procedure_Redd(env_b,D).Save_matrix(Data_Way + base_file_name + "_env_bot_redd.txt");//cout << "Redd\tb0\n";
        ///bot sm3
        procedure_Redd(env_b3,D).Save_matrix(Data_Way + base_file_name + "_env_bot_sm3_redd.txt");//cout << "Redd\tb3\n";
        ///bot sm10
        procedure_Redd(env_b10,D).Save_matrix(Data_Way + base_file_name + "_env_bot_sm10_redd.txt");//cout << "Redd\tb10\n";
        ///topbot sm0
        procedure_Redd(env_tb,D).Save_matrix(Data_Way + base_file_name + "_env_tb_redd.txt");
        ///topbot sm3
        procedure_Redd(env_t3b3,D).Save_matrix(Data_Way + base_file_name + "_env_tb_sm3_redd.txt");
        ///topbot sm10
        procedure_Redd(env_t10b10,D).Save_matrix(Data_Way + base_file_name + "_env_tb_sm10_redd.txt");

 /**
    ofstream fout_ampl( Data_Way + base_file_name + "_ampl_res.txt");
        ///анализ амплитуды
        ///найти точку перегиба огибающей
        matrix d_env_top;///производная огибающей
            d_env_top = derevative_3_dot_h(env_t, 1,0);///ПРОИЗВОДНАЯ БЕРЕТСЯ В ПРЕДПОЛОЖЕНИИ ПОСТОЯННОГО ШАГА
                d_env_top.abs_this_column(1); ///модуль столбца, чтобы искать экстремум
        int index_mid_angle = d_env_top.get_max_elind_column(1);///точка экстремума - точка перегиба

        d_env_top.column_multiply(1,1.0/(d_env_top.get_data(index_mid_angle,1))); ///нормируем на максимум
            d_env_top.Save_matrix(Data_Way + base_file_name + "_env_top_derev.txt");

        double mid_angle = env_t10.get_data(index_mid_angle, 1); ///угол точки перегиба
                    cout << "index_mid_angle = " << index_mid_angle <<"\tmid_angle = " << mid_angle << endl;

        ///поиск полки
        int p_back, p_fow;
        double plateau_condition = 0.1;
        cout << "PLATEAU 1\t" << plateau_condition << endl;
        p_back = get_plateau_backward(env_t10, 0, 1, plateau_condition);
        p_fow = get_plateau_forward(env_t10, 0, 1, plateau_condition);
            cout << "p_back = " << p_back << "\tp_fow = " << p_fow << endl;
            cout << env_t10.get_data(p_back, 1) << "\t" << env_t10.get_data(p_fow, 1);

        cout << "PLATEAU 2\n";
        ///сместили на мид англ (середина, где сменяется характер роста, в нуле)
        env_t10.column_shift(1,-1.0*mid_angle);
        ///ищем полки
        p_back = get_plateau_backward(env_t10, 0, 1, 0.1);
        p_fow = get_plateau_forward(env_t10, 0, 1, 0.1);
            cout << "p_back = " << p_back << "\tp_fow = " << p_fow << endl;
            cout << env_t10.get_data(p_back, 1)+mid_angle << "\t" << env_t10.get_data(p_fow, 1) + mid_angle;///вернем смещение

    fout_ampl.close();*/
};

void procedure_ampl_analisys(string Data_Way, string base_file_name, oscillation& D, double freq, double angle_shift)
    {procedure_ampl_analisys(Data_Way, base_file_name, "", D, freq, angle_shift);}

matrix data_cutter(matrix A, int column_num, double column_data1, double column_data2)
{cout << "duta_cutter" << endl;//A.info();
    matrix res;
    int index_start = 0, index_end = 0;//cin >> index_end;
        //cout << column_data1 << "\t" << column_data2 << endl;
        while(A.get_data(index_start, column_num)<=column_data1)
                index_start++;
        index_end = index_start;//cout << A.get_data(index_end, 0) << "\t" << column_data2 << "\t" << index_end << endl;
        if(A.get_data(A.get_height()-1,column_num) >= column_data2)
        {//cout << A.get_data(index_end, 0) << "\t" << column_data2 << "\t" << index_end << endl;
            while(A.get_data(index_end, column_num)<=column_data2)
                index_end++;
        }else{
            index_end = A.get_height()-1;
        }       //cout << A.get_data(index_end, 0) << "\t" << column_data2 << "\t" << index_end << endl;
    res = A.get_matrix_part(index_start,index_end,0,A.get_width());res.info();
    return res;
}

void Redd_Cut_file(string Data_Way, string instruct_name)
{
    string new_name, out_name = "";
    matrix A;
    double time_end = 0, time_start = 0;

    ifstream fin(Data_Way + instruct_name + ".txt");
    while(new_name != "START")
        fin >> new_name;///buff
            fin >> new_name; ///input file_name

    while(new_name != "END")
    {
        fin >> time_start >> time_end >> out_name; cout << "Part = " << time_start << " <-> " << time_end << endl;

        if(FileExists(Data_Way + new_name + "_env_top_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_top_redd.txt");//A.info();cout<< A;
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_top_0.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_bot_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_bot_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_bot_0.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_top_sm3_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_top_sm3_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_top_3.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_bot_sm3_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_bot_sm3_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_bot_3.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_top_sm10_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_top_sm10_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_top_10.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_bot_sm10_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_bot_sm10_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_bot_10.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_tb_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_tb_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_tb_0.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_tb_sm3_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_tb_sm3_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_tb_3.txt");
        }

        if(FileExists(Data_Way + new_name + "_env_tb_sm10_redd.txt"))
        {
            A.Load_matrix(Data_Way + new_name + "_env_tb_sm10_redd.txt");
                data_cutter(A, 0, time_start, time_end).Save_matrix(Data_Way + out_name + "_Redd_tb_10.txt");
        }

    fin >> new_name;
    }
    fin.close();
}



double inertia_cg(double xl, double yl, double h, double m1, double m2, double alpha1, double alpha3)
{
    double r2 = sqrt(pow((xl+h),2) + pow(yl, 2));
        return r2*m2/m1*sin(alpha3*0.0174533)/sin(alpha1*0.0174533);
}

double inertia_cg1(double xl, double yl, double h, double m1, double m2, double a_free, double a_aft, double a_w)
{
    double r2 = sqrt(pow((xl+h),2) + pow(yl, 2));
    double gamma, a0, a1, a3;
    gamma = asin(yl/r2)/0.0174533;
    a1 = fabs(a_free - a_w);
    a0 = fabs(a_free - a_aft); //a0 = 3.81;
    a3 = 90 - gamma - a0 - a1;
    cout << "r2 = " << r2 << "\t" << m1 << "\t" << m2 << "\t" << gamma << "\t" << a0 << "\t" << a1 << "\t" << a3 << "\t";
    return r2*m2/m1*sin(a3*0.0174533)/sin(a1*0.0174533);
}

double inertia_momentum(double xl, double yl, double h, double m1, double m2, double a_free, double a_aft, double a_w, double g, double T)
{///ЗДЕСЬ НЕ УЧИТЫВАЕТСЯ m1 !!!
    //inertia_cg1(xl, yl, h, m1, m2, a_free, a_aft, a_w);
    double r2 = sqrt(pow((xl+h),2) + pow(yl, 2));
    double gamma, a0, a1, a3;
        gamma = asin(yl/r2)/0.0174533;
        a1 = fabs(a_free - a_w);
        a0 = fabs(a_free - a_aft);//a0 = 3.76; //a0 = 72.09; //a0 = 3.76;//a0 = 3.81; по старым измерениям //a0 = 3.76; подобранный
        a3 = 90 - gamma - a0 - a1;
        ///cout << gamma << "\t" << a0 << "\t" << a1 << "\t" << a3 << "\t";
        //a3 = 90 - gamma - a0 - a1;
        m2 *=0.001;
        r2 *= 0.001;
            return g*r2*m2*sin(a3*0.0174533)/sin(a1*0.0174533)*pow(T,2)/4/pow(PI,2);
}

void procedure_inertia_analysis(double Isw, double m1, double m2, double h, double d, double alpha0, double alpha1, double alpha3, double yl, double r2, double g, double T)
{
    double gamma, r1, r3, I1, Is, Is_exp;
    ///нужно для дальнейшего счета
    r1 = r2*m2/m1*sin(alpha3*0.0174533)/sin(alpha1*0.0174533);
        cout<< "r1 = " << r1 << "\tr2 = " << r2 << "\tm2 = " << m2 << "\tm1 = " << m1 << "\talph3 = " << alpha3 << "\talpha1 = " << alpha1 << endl;
    gamma = asin(yl/r2/1000.0)*180.0/PI;

    ///проверка
    r3 = (r1*sin(alpha0*0.0174533)*m1 + r2*cos(gamma*0.0174533)*m2)/(m1+m2)/sin((alpha0 + alpha1)*0.0174533);///через X
        cout << "r3 = " << r3 << endl;

    I1 = (3.0*d/2*d/2+h*h)*m2/12.0;///циллиндр
    Is = Isw + I1 + m2*r2*r2;
    Is_exp = (m1+m2)*g*r3*pow(T,2)/4/pow(PI,2);
        cout <<"Is = " << Is << "\tIsw = " << Isw << "\tI1 = " << I1 << "\tmr2 = " << m2*r2*r2 << endl;
            cout << "Is_exp = " << Is_exp << "\tIs/Is_exp = " << 1 - Is/Is_exp << endl;
}

void procedure_Mz_FFT(string Data_Way, string base_file_name, double time_from, double time_to, double window_step_time, int window_size_power, int window_zeroes_power )
{cout << "procedure_Mz_FFT" << endl;
    matrix m; m.Load_matrix(Data_Way + base_file_name + ".txt");///загрузка исходных данных
        matrix Data, Data_part, fft,
            time_freq,///результирующие данные по главной гармонике
                t_f_buff(1,2);///частота и время в окне
    Data = m.get_column(0);///time
        Data = Data.merge_width(m.get_column(1), Data.get_width());///angle
    ///служебные
    int window_size = pow(2, window_size_power);///размер данных окна ///мб делать любого размера?
            int fft_num = 0;///счетчик окна
                int window_step_num = window_step_time/(m.get_data(1,0) - m.get_data(0,0)); ///шаг смещения окна
                    cout << "window_size = " << window_size << "\twindow_step_num = " << window_step_num << "\ttotal_length = " << m.get_height() << endl;

    int h1 = 0, h2 = h1 + window_step_num;///modify time_from, time_to ///начальная и конечная точка для обработки

    while(h2 < Data.get_height()) ///modify time_to
    {//cout << "h1 = " << h1 << "\th2 = " << h2 << endl;
        Data_part = Data.get_matrix_part(h1, h2, 0, 2);//cout << Data_part << endl;
            h1 += window_step_num;
                h2 = h1 + window_step_num;

        fft = Get_Spectr_column(Data_part, 0, 1, window_zeroes_power);///fft_height = 2^(window_size_power + window_zeroes_power)
            //cout << "SPECTR\n";
        t_f_buff.set_data(0, 0, Data_part.get_data(0,0));///time//cout << t_f_buff << endl;
            t_f_buff.set_data(0, 1, fft.get_data(get_max_element_index_from_matrix_column(fft, 1, 0, (int)(fft.get_height()/2)), 0));///freq//cout << "max_freq = " << t_f_buff << endl;
                time_freq = time_freq.merge_height(t_f_buff, time_freq.get_height());

        fft.Save_matrix(Data_Way + base_file_name + "_fft_" + manager_int_num(fft_num) + ".txt");
            fft_num++;
    }
    time_freq.Save_matrix(Data_Way + base_file_name + "_main_freq" + ".txt");
}

/// DIFFERENCIAL EQUATIONS

///IMPLICIT
/*
void Calc_VdP_equation_4(string Data_Way, string base_file_name, int fix_max)
{
    int N_steps;
        ///ВИД УРАВНЕНИЯ
        double T, h, lam, mu;
        double gb, gc, pb, pc, x0, y0, z0, shift_angle = 0, shift_time = 0;
        matrix C_RES, RES, d_RES, Env_RES, Ampl_RES(1,4), Empt;///Ampl_RES //lam, mu, rho-, rho+

        string s_buff = "", new_name;
        ifstream fin(Data_Way + "VdP_4_param.txt");

    while(new_name != "START")
        fin >> new_name;///buff
    fin >> new_name; ///data
    cout << "new_name = " << new_name << endl;
    while(new_name != "END")
    {
            fin >>
                T >> h >>
                    x0 >> y0 >> z0 >> ///начальные значения в градусах!
                            lam >> mu;/// коэффициенты уравнения при производной порядка 1, 1 и (x*x), 0 (при 2 = 1)
            N_steps = T/h;///количество шагов по времени и указанному шагу

        ///d2y/dx2-(lam+mu*y2-y4)dy/dx+x
        ///d2y/dx2+(-lam-mu*y2+y4)dy/dx+x
        pb = 1.0; pc = 0.0;
        C_RES = equation_func_VdP_4(Data_Way, base_file_name, h, N_steps, fix_max, -lam, -mu, 1.0, 1, pb, pc, x0, y0, z0);;
            C_RES.Save_matrix(Data_Way + new_name + "_VdP_eq.txt");

        RES = C_RES.get_column(0);
            RES = RES.merge_width(C_RES.get_column(1), RES.get_width());
        d_RES = C_RES.get_column(0);
            d_RES = d_RES.merge_width(C_RES.get_column(2), RES.get_width());

        Env_RES = Get_Envelope_top(RES, d_RES);
            Env_RES.Save_matrix(Data_Way + new_name + "_envelop_eq.txt");
        ///амплитуда минимальная и максимальная
        Ampl_RES.set_data(0,0, lam);
        Ampl_RES.set_data(0,1, mu);
        Ampl_RES.set_data(0,2, Env_RES.get_min_el_column(1));
        Ampl_RES.set_data(0,3, Env_RES.get_max_el_column(1));
            Ampl_RES.Save_matrix_end(Data_Way + new_name + "_ampl_eq.txt");

        fin >> new_name ;/// следующий или конец?
        C_RES = Empt;
    }
    fin.close();
}

void Calc_VdP_equation_4_mdyn(string Data_Way, string base_file_name, string mode, int fix_max)
{cout << "Calc_VdP_equation_4_mdyn(string Data_Way, string base_file_name, string mode, int fix_max)\n";
    int N_steps;
        ///ВИД УРАВНЕНИЯ
        double T, h, e0, e2, e4, mdyn, gc0, gc2, gc4; ///время счета, шаг счета, коэффициенты.
        double gb, gc, gd, pb, pc, pd, x0, y0, z0, shift_angle = 0, shift_time = 0;
        matrix C_RES, RES, d_RES, Env_RES, Ampl_RES(1,4), const_column;///Ampl_RES //lam, mu, rho-, rho+

        string s_buff = "", new_name;
        ifstream fin(Data_Way + "VdP_4_mdyn_param.txt");

    while(new_name != "START")
        fin >> new_name;///buff
    fin >> new_name; ///data
        cout << "DATA = " << new_name << endl;
    while(new_name != "END")
    {
            fin >>  T >> h >>
                    x0 >> y0 >> z0 >> ///начальные значения в градусах!
                            e0 >> e2 >> e4 >> mdyn >>
                                gc0 >> gc2 >> gc4;/// коэффициенты уравнения при производной порядка 1, 1 и (x*x), 0 (при 2 = 1)
            N_steps = T/h;///количество шагов по времени и указанному шагу

        gb = 0.0; gc = -1.0; gd = 0.0; ///dz/dx = gb()*z - x
        pb = 1.0; pc = 0.0; pd = 0.0; ///dy/dx = z

        if (mode == "0")
            C_RES = equation_func_VdP_4_mdyn(Data_Way, base_file_name, h, N_steps, fix_max, e0, 1.0*e2, 1.0*e4, gc, gd, pb, pc, pd, mdyn, x0, y0, z0);
        if (mode == "exp_eq")
            C_RES = equation_func_VdP_4_mdyn(Data_Way, base_file_name, h, N_steps, fix_max, e0, 1.0*e2, 1.0*e4, gc0, gd, pb, pc, pd, mdyn, x0, y0, z0);
        if (mode == "VdP4_mdyn_poly_mst")
            C_RES = equation_func_VdP_4_mdyn_poly_mst(Data_Way, base_file_name, h, N_steps, fix_max,e0, 3.0*e2, 5.0*e4,gc0, gc2, gc4, gd, pb, pc, pd, mdyn,x0, y0, z0);
        if (mode == "advanced_RK4")
        {
            cout << "advanced_RK4\n\tf = Poly[4], g = Poly[1]\n";
                C_RES = equation_func_polynomial4order_f_g(Data_Way, base_file_name, h, N_steps, fix_max,
                                                   [](double fa0, double fa2, double fa4, double f_arg){return fa0 + fa2*pow(f_arg,2) + fa4*pow(f_arg,4);},
                                                   [](double ga0, double ga2, double ga4, double g_arg){return ga0 + ga2*pow(g_arg,2) + ga4*pow(g_arg,4);},
                                                   e0, 3.0*e2, 5.0*e4,gc0, gc2, gc4,
                                                   gd, pb, pc, pd, mdyn, x0, y0, z0);
        }

        ///сохранение результатов численного решения уравнения и сч матрицы Якоби: x,y,z,lam1Re,lam1Im,lam2Re,lam2Im
        C_RES.Save_matrix(Data_Way + new_name + "_VdP_eq" + ".txt");

        ///сохранение в форме oscillation данных для обработки: count1, flow3, angle1, dangle1, freq1, time1 = 8
        const_column.make_matrix_data(C_RES.get_height(),1,1.0);/// имеется время, угол, скорость
        RES = RES.merge_width(const_column);///count
        RES = RES.merge_width(const_column);///flow1
        RES = RES.merge_width(const_column);///flow2
        RES = RES.merge_width(const_column);///flow3
        RES = RES.merge_width(C_RES.get_column(1));///angle
        RES = RES.merge_width(C_RES.get_column(2));///dangle
        RES = RES.merge_width(const_column);///freq
        RES = RES.merge_width(C_RES.get_column(0));///time
        RES.Save_matrix(Data_Way + new_name + "_r" + ".txt");

        C_RES.delete_data();///освобождение для следующей задачи
        fin >> new_name ;/// следующий или конец?
    }
    fin.close();
}*/

void Calc_VdP4_FrontRear(string Data_Way, string base_file_name, string mode, int fix_max)
{
    ///данные
    matrix C_RES, RES, d_RES, Env_RES, Ampl_RES(1,4), Empt;

    ///*********************
    ///счет уравениня как совокупности коэффициентов для донной фрональной и донной частей модели
    double e0, e2, e4; /// параметры донной части
        double time_lag;    ///временная задержка (коэффициент при функции mdyn)
    double Mst_front; ///коэффициент при фронтальной части
    double gc0, gc2, gc4, gb, gc, gd, pb, pc, pd, x0, y0, z0; ///параметры системы уравнений
         double T, h;
            int N_steps;
    ///**********************


    string s_buff = "", new_name; ///буфер считывания и имя решения
    ifstream fin(Data_Way + "VdP4_FrontRear.txt");

        while(new_name != "START")///исполянются строки начиная со START
            fin >> new_name;///buff
        fin >> new_name; ///data ///имя решения
            cout << new_name << endl;
        while(new_name != "END")/// исполнение закаечивается при начале строки со STOP
        {

            ///READ INIT
            ///загрузка данных из файла
            fin >>
                T >> h >> ///время для рассчета и шаг по времени
                    x0 >> y0 >> z0 >> ///начальные условия x = t, y = alpha, z = d(alpha)/dt
                            e0 >> e2 >> e4 >> time_lag >> ///time_lag безразмерный
                                 Mst_front;
            N_steps = T/h;///количество шагов по времени и указанному шагу
            ///READ INIT ///

        ///разложение для статической части (предполагается что передняя часть имеет лишь линейную составляющую момента!)
            gc0 = e0 + Mst_front;
            gc2 = e2;
            gc4 = e4;

        ///решение сведением к СЛУ
        ///dz/dx = gb*z+gc*y+gd*x
        ///dy/dx = pb*z+pc*y+pd*x
        ///gb фактически это e0, e2, e4
        gb = 0.0; gc = gc0; gd = 0.0;
        pb = 1.0; pc = 0.0; pd = 0.0;
            C_RES = equation_func_VdP_4_mdyn_poly_mst(Data_Way, base_file_name, h, N_steps, fix_max, e0, e2, e4, gc0, gc2, gc4, gd, pb, pc, pd, time_lag, x0, y0, z0);

            RES = C_RES.get_column(0);
                RES = RES.merge_width(C_RES.get_column(1), RES.get_width());
            d_RES = C_RES.get_column(0);
                d_RES = d_RES.merge_width(C_RES.get_column(2), RES.get_width());
        ///постобработка
            ///построение огибающей
            Env_RES = Get_Envelope_top(RES, d_RES);
                Env_RES.Save_matrix(Data_Way + new_name + "_envelop_eq.txt");

        ///ШАГ итерации
        fin >> new_name ;/// следующий или конец?
            ///очистка данных
            C_RES = Empt;
        }
    fin.close();
}


void Analys_VdP4_mdyn(string Data_Way, string base_file_name, string name_add, string angle_val, double mdyn, double epsilon4, double epsilon2 )
{cout << "Analys_VdP4_mdyn()\n";
    ///epsilon2 = 3.0*(A1-A2)*A2/A1;
    ///epsilon4  = 5.0*A2/A1;
    cout << "Analys_VdP4_mdyn\n = " << Data_Way + base_file_name << endl;
    cout << "dy2/dy2 + " << mdyn << "*("<<epsilon4 << "*x^4 + " << epsilon2 << "*x^2 + " << 1 << ") - x" << endl;

    matrix buff;
        buff.Load_matrix(Data_Way + base_file_name + ".txt");
        cout << "BUFF INFO\n"; buff.info();
        cout << "time/coordinate/...etc...\n";
    matrix data;    ///t //x // dx/dt
        data = buff.get_column(0);//data.info();
            data = data.merge_width(buff.get_column(1), data.get_width());

    for(int i = 0; i < data.get_height(); i++)
    {
        data.set_data(i,0, data.get_data(i,0));
        if (angle_val == "deg")
            data.set_data(i,1, data.get_data(i,1)/180.0*PI);
    }

    matrix d1,d2;
        d1 = derevative_3_dot(data,1,0);
            d2 = derevative_3_dot(d1,1,0);
    d1.erase_column(0); d2.erase_column(0);
    ofstream fout(Data_Way + base_file_name + "_an_md.txt");
    //matrix mdyn(data.get_height(), 2);
        for(int i = 0; i < data.get_height(); i++)
        {
            //fout << i << "\t" << mst*v/l/w0 << "\t" << (d2.get_data(i,0)-data.get_data(i,1))/(epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0)/d1.get_data(i,0) << endl;
          //  mdyn.set_data(i, 0, mst*v/l/w0*(d2.get_data(i,0)-data.get_data(i,1))/(epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0)/d1.get_data(i,0));
            //cout << d2.get_data(i,0) << "\t" << data.get_data(i,1) << "\t" << d1.get_data(i,0) << "\t" << 0.01*(epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0) << endl;
          //  mdyn.set_data(i, 0, d2.get_data(i,0)-data.get_data(i,1)-0.01*(epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0)*d1.get_data(i,0));
            //cout << data.get_data(i,0) << "\t" << data.get_data(i,1) << "\t" << d1.get_data(i,0) << "\t" << d2.get_data(i,0) << (epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0) << "\t" << (epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0)*d1.get_data(i,0) << endl;

            fout << data.get_data(i,0) << "\t" << data.get_data(i,1) << "\t" << d1.get_data(i,0) << "\t" <<  d2.get_data(i,0) << "\t" <<  (epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0) << "\t"
                << "\t" << d2.get_data(i,0) - data.get_data(i,1)  << "\t" << + mdyn*(epsilon4/mdyn*pow(data.get_data(i,1),4) + epsilon2/mdyn*pow(data.get_data(i,1),2) + 1.0)*d1.get_data(i,0)
                    << "\t" << d2.get_data(i,0) - data.get_data(i,1) * mdyn/(epsilon4*pow(data.get_data(i,1),4) + epsilon2*pow(data.get_data(i,1),2) + 1.0)/d1.get_data(i,0) << endl;
         }
    fout.close();
            //mdyn.set_data(i, 0, w0*sqrt(-1.0*iz*mst)*(d2.get_data(i,0)-w0*w0*data.get_data(i,1))/(1.0+(A1+A2)*pow(data.get_data(i,1),2)+A1*A2*pow(data.get_data(i,1),4))/d1.get_data(i,0));

    data = data.merge_width(d1, data.get_width());
    data = data.merge_width(d2, data.get_width());
    data.Save_matrix(Data_Way + "data.txt");
   // mdyn.Save_matrix(Data_Way + base_file_name + name_add + "_mdyn.txt");
}

/*void Calc_Equation_procedure()
{

}*/

void Get_Cycle_Radii(string Data_Way, string file_name, double from, double to, double step)
{cout << "Get_Cycle_Radii\n" << endl;
    int N = 0;
    double lam = from, mu = from;
    N = (to - from)/step;
        N++;///переход через ноль
            cout << "N = " << N << endl;
    matrix R(N*N+N,4);///mu lam rho+ rho- ///+N - один на конечный симол для каждого куска
    R.info();
    for(int i = 0; i < N*N; i++)
    {//cout << "i = " << i << "\t";
        for(int k = 0; k < N+1; k++ )
        {//cout << "k = " << k << "\t";
            R.set_data(i, 0, mu);
            R.set_data(i, 1, lam);
            if((mu*mu+8.0*lam)>=0)
                {
                    R.set_data(i, 2, (mu+sqrt(mu*mu+8.0*lam))/4.0);
                    R.set_data(i, 3, (mu-sqrt(mu*mu+8.0*lam))/4.0);
                }else{
                    R.set_data(i, 2, 0);
                    R.set_data(i, 3, 0);
                }
            i++;
            lam +=step;
        }lam = from;
        //cout << endl;
            mu +=step;
    }
    //R.erase_column(3);
    R.Save_matrix(Data_Way + file_name + ".txt");
};

///FLUENT REPORT
void get_mz_results(string Data_Way, string base_name, string result_name, int from_index, int to_index)
{
    int rh = 0, fh = 0, r_counter = 0;
    ifstream fin;
    string s_buff = "";
    double d_buff = 0;
    ///пересчитать количество существующих файлов
    for(int i = from_index; i <= to_index; i++)
        {
            s_buff = Data_Way + "/" + manager_int_num(i) + "degree/" + base_name;
                if(FileExists(s_buff.c_str()))///no format
                    rh++;
        }
    ///создали нужного размера
    matrix result(rh, 2); /// angle , mz убрано iteration

    ///по каждому файлу
    for(int angle_num = from_index; angle_num <= to_index; angle_num++) ///по карте читать чтобы дважды не вызывать FileExists
        ///чтение последней строки (double,double, dobuble) из некоторого файла файла
        if(FileExists((Data_Way + "/" + manager_int_num(angle_num) + "degree/" + base_name).c_str()))///no format
        {
            ///открыли нужный файл
            fin.open((Data_Way + "/" + manager_int_num(angle_num) + "degree/" + base_name).c_str()/*, ios_base::openmode mode = ios_base::in*/);
            ///подсчет полного количества строк в файле
            fh = 0;
            while(!fin.eof())
            {
                getline(fin,s_buff);
                    fh++;
            }
            //fin.close();
            ///вернулись в начало
            fin.seekg (0, fin.beg);
            ///получение нужных данных из файла
            //fin.open(Data_Way + "/" + manager_int_num(angle_num) + "degree/" + base_name/*, ios_base::openmode mode = ios_base::in*/);
            for(int j = 0; j < fh; j++)
            {
                if(j == fh-1-3)///последняя пустая, предпоследняя нам нужна
                {
                    //cout << "j = \t" << j << "\t" << "j = fh - 1\t" << j << "\t" << s_buff << endl;
                    result.set_data(r_counter, 0, angle_num);
                    fin >> d_buff;///итерация
                        //result.set_data(r_counter, 1, d_buff);//cout << d_buff << "\t";
                    fin >> d_buff;///значение
                        result.set_data(r_counter, 1, d_buff);//cout << d_buff << "\n";
                    j = fh;
                    r_counter++;
                }
                getline(fin,s_buff);///читаем все не нужные строки
            }
            fin.close();
        }
    result.Save_matrix(Data_Way + result_name +".txt");
}

void get_mz_results1(string Data_Way, string base_name, string result_name, int from_index, int to_index)
{
    int rh = 0, fh = 0, r_counter = 0;
    ifstream fin;
    string s_buff = "";
    double d_buff = 0;
    ///пересчитать количество посчитанных
    for(int i = from_index; i <= to_index; i++)
    {
        s_buff = Data_Way + "/" + manager_int_num(i) + "degree/" + base_name;
            if(FileExists(s_buff.c_str()))///no format
                rh++;
    }
    ///создали нужного размера
    matrix result(rh, 3); /// angle, iteration, mz
    ///по каждому файлу
    for(int angle_num = from_index; angle_num <= to_index; angle_num++)
        ///чтение последней строки (double,double, dobuble) из некоторого файла файла
        if(FileExists((Data_Way + "/" + manager_int_num(angle_num) + "degree/" + base_name).c_str()))///no format
        {
            fin.open((Data_Way + "/" + manager_int_num(angle_num) + "degree/" + base_name).c_str());
                ///последняя
                fin.seekg (0, fin.end);
                ///количество позиций!!!
                fh = fin.tellg();
                cout << "height = " << fh << endl;
                ///предпоследняя
                fin.seekg (0, fin.beg);
                ///читаем
                    result.set_data(r_counter, 0, angle_num);
                fin >> d_buff;///итерация
                    result.set_data(r_counter, 1, d_buff); cout << d_buff << "\t";
                fin >> d_buff;///значение
                    result.set_data(r_counter, 2, d_buff); cout << d_buff << "\n";
                r_counter++;
            fin.close();
        }
    result.Save_matrix(Data_Way + result_name +".txt");
}
