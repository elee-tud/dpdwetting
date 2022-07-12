#ifndef __TIMER__HPP
#define __TIMER__HPP
#include "types.hpp"
#include <iostream>
#include <iomanip>
#include <ctime>
#include <string>
#include <sstream>
namespace dpd{
class Timer{
private:
	time_t start_time;
	int total_steps;
    int start_step;
	time_t program_time;
	time_t setu_time;
	time_t read_time;
	time_t calc_time;
	time_t comm_time;
	time_t program_start;
	time_t readt_start;
	time_t calct_start;
	time_t commt_start;
	time_t setut_start;
	int mpirank;
    std::string remaintime;

public:
	Timer(){};
	Timer(int _mpirank, int _total_steps);
	~Timer();
    void reset(int _total_steps);
    void reset(int _start_step, int _final_steps);
	void printProgressRemainingTime(int cur_step);
	void printProgressRemainingTime(int cur_step, real max_force);

	void printFinalTime();
    std::string getRemainingTime(){return remaintime;}


};

struct Time{
	long hrs;
	long mins;
	long secs;
};

Time convertSecondsToTime(time_t seconds);
std::string convertTimeToString(Time t);
std::string convertSecondsToString(time_t seconds);

};

#endif

