import sys, datetime, time
import numpy as np

BLOCK = ['', chr(0x258F), chr(0x258E), chr(0x258D), chr(0x258C), chr(0x258B), chr(0x258A), chr(0x2589), chr(0x2588)]

def block_string(_per, _bar_length=20):
    num1 = int(_per * _bar_length)
    remain = _per * _bar_length - num1
    remain = max(int(remain * 8), 0)

    blk = f"{chr(0x2588) * num1 + BLOCK[remain]:<{_bar_length}}"

    return blk

class Status_Bar():
    def __init__(self, _times=None):
        self.total_itr = _times
        self.time_sum = 0
        self.present = 0
        self.prop_name = "None"
        self.time_len = len(str(self.total_itr))

        self.verbose = False

        self._bar_length = 20

        try:
            self.total_itr = int(self.total_itr)
        except:
            pass

    def set_total_times(self, _times):
        self.total_itr = int(_times)
        self.time_len = len(str(self.total_itr))

    def add_total_times(self, _times):
        self.total_itr += _times

    def sub_total_times(self, _times):
        self.total_itr -= _times

    def set_prop_name(self, prop_name):
        self.prop_name = prop_name

    def set_present(self, _times):
        self.present = _times

    def _predict_finish_time(self):
        if self.present == 0:
            return datetime.datetime(2999, 12, 31, 23, 59, 59)
        current_time = datetime.datetime.now()
        avg_time = self.time_sum / self.present
        remain_time = datetime.timedelta(seconds=avg_time * (self.total_itr - self.present))

        return current_time + remain_time

    def output(self):
        if self.verbose:
            sys.stderr.write('\r')

            current_per = self.present / self.total_itr
            output_str = f"{self.prop_name} [{block_string(current_per, self._bar_length)}] {self.present:>{self.time_len}}/{self.total_itr:0>{self.time_len}} ({current_per:.3%}), ETC: {self._predict_finish_time():%Y/%m/%d,%H:%M:%S}  "

            sys.stderr.write(output_str)
            # sys.stdout.flush()

    def append(self, _time):
        self.present += 1
        self.time_sum += _time

    def finish(self):
        sys.stderr.write('\n')