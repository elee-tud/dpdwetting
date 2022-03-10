#include "timer.hpp"
using namespace dpd;
Timer::Timer(int _mpirank, int _total_steps):mpirank(_mpirank), total_steps(_total_steps){
	if(mpirank==0){
		program_start=time(NULL);
        start_step=0;
	}

}
Timer::~Timer(){}

void Timer::reset(int _total_steps){
    if(mpirank==0){
        total_steps=_total_steps;
        program_start=time(NULL);
    }
    return;
}

void Timer::reset(int _start_step, int _final_step){
    if(mpirank==0){
        start_step=_start_step;
        total_steps=_final_step;
        program_start=time(NULL);
    }
    return;
}

void Timer::printProgressRemainingTime(int cur_step){
	if(mpirank==0 && cur_step%10==0){
		real done=static_cast<real>(cur_step-start_step)/(total_steps-start_step);
		time_t cur_time=time(NULL);
		time_t elapsed=cur_time-program_start;
		time_t remain=(elapsed/done*(1.-done));
        remaintime=convertSecondsToString(remain);
		

        std::cout << "\r Step "<< std::setw(10) << cur_step <<"(" << std::setw(4) << std::right << std::setprecision(3) << std::right <<  done*100 <<"%) out of " << total_steps << " is done. Remaining time is " << remaintime << std::flush;
        if(cur_step==total_steps)
            std::cout << std::endl;
	}
}

void Timer::printProgressRemainingTime(int cur_step, real max_force){
	if(mpirank==0 && cur_step%1==0){
		real done=static_cast<real>(cur_step)/total_steps;
		time_t cur_time=time(NULL);
		time_t elapsed=cur_time-program_start;
		time_t remain=(elapsed/done*(1.-done));
	

        std::cout << "\r Step "<< std::setw(10) << cur_step <<"(" << std::setw(4) << std::right << std::setprecision(3) << std::right <<  done*100 <<"%) out of " << total_steps << " is done. Current maxium force = " << max_force << ". Remaining time is " << convertSecondsToString(remain) << std::flush;
        if(cur_step==total_steps)
            std::cout << std::endl;
    }
}


void Timer::printFinalTime(){
	if(mpirank==0){
        program_time=time(NULL)-program_start;
        std::cout << std::endl;
              std::cout << "***************************************************************************"<< std::endl;
		      std::cout << std::setw(30) << "Total calc. Time"<< "=" << convertSecondsToString(program_time) << std::endl;
		      std::cout << std::setw(30) << "Time/step"<< "=" << convertSecondsToString(static_cast<time_t>(program_time/(total_steps-start_step))) << std::endl;
              std::cout << "***************************************************************************"<< std::endl;
	}
}


Time dpd::convertSecondsToTime(time_t seconds){
	Time t;
	t.hrs=static_cast<long>(seconds)/3600;
	t.secs=static_cast<long>(seconds)%3600;
	t.mins=t.secs/60;
	t.secs=t.secs%60;
	return t;
}

std::string dpd::convertTimeToString(Time t){
    std::stringstream out;
    out << std::setw(5) << std::right << t.hrs << " hrs " 
	<< std::setw(2) << std::right << t.mins << " min "
	<< std::setw(2) << std::right << t.secs << " sec.";
	return out.str();
}

std::string dpd::convertSecondsToString(time_t seconds){
	return convertTimeToString(convertSecondsToTime(seconds));
}
	
