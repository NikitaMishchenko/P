#ifndef MY_TIMER_H_INCLUDED
#define MY_TIMER_H_INCLUDED

#include <thread>
#include <chrono>


void time_func(int sleep_ms)
{std::cout << "time_func(int)\n";
    //this_thread::sleep_for(chrono::milliseconds(sleep_ms));
}

#endif // MY_TIMER_H_INCLUDED
