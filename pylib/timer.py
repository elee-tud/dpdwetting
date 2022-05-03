import sys
import time
import datetime

class Timer:
    def __init__(self, _total_steps):
        self.total_steps=_total_steps
        self.start_time=time.time()
    def process_done(self, curr_step):
        curr_step += 1
        curr_time=time.time()
        used_time=curr_time-self.start_time
        estimated_finish=curr_time+used_time*float(self.total_steps-curr_step)/float(curr_step)
        print('Current step %d out of total steps %d --> Estimated time to finish: %s'%(curr_step, self.total_steps, time.asctime(time.localtime(estimated_finish))), end="\r", flush=True)
    def finished(self):
        curr_time=time.time()
        print("\nCalculation finished")
        print("Time consumed: %s (%s per step)\n"%(time.strftime("%H hrs %M min %S sec", time.gmtime(curr_time-self.start_time)), time.strftime("%H hrs %M min %S sec", time.gmtime((curr_time-self.start_time)/self.total_steps))))


