/*
 * cross_timer.h
 *
 *  Created on: Apr 24, 2017
 *      Author: shurbevski
 */
// This is the file used to measure running time.

#ifndef CROSS_TIMER_H_
#define CROSS_TIMER_H_

#include <chrono>

// Chrono based time measurement
class Timer {
public:
    Timer() {
        reset();
    }
    void reset() {
        m_timestamp = std::chrono::high_resolution_clock::now();
    }
    double diff() {
        std::chrono::duration<double>
        	fs = std::chrono::high_resolution_clock::now()
        		- m_timestamp;
        return fs.count();
    }
private:
    std::chrono::high_resolution_clock::time_point m_timestamp;
};



struct globalTimeKeeper {
	static Timer tt;
	static bool enabled;
	static float timeOut;
};
#endif /* CROSS_TIMER_H_ */
