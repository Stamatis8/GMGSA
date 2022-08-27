#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iostream>

struct timer {

	timer() {};

	void begin() {
		this->start_time = std::chrono::high_resolution_clock::now();
	}

	void display() {
		std::cout << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - this->start_time).count();
	}

private:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};//timer

#endif// TIMER_HPP
