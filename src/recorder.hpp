#pragma once

#include <string>
#include <map>
#include <chrono>
#include <iostream>

class Recorder {
public:
	class Timer {
		std::chrono::time_point<std::chrono::steady_clock> start;
		std::string name;

	public:
		Timer(std::string name) : name(name), start(std::chrono::steady_clock::now()) {}
		~Timer() {
			Recorder::submit(name, std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - start));
		}
	};

private:
	inline static std::map<std::string, std::chrono::microseconds> timers;
	inline static std::map<std::string, int> calls;

public:
	inline Recorder() = delete;
	inline ~Recorder() = delete;

	static inline void  submit(std::string name, std::chrono::microseconds time) {
		if (timers.find(name) == timers.end()) {
			timers.emplace(name, time);
			calls.emplace(name, 1);
		}
		else {
			timers[name] += time;
			calls[name]++;
		}
	}

	static inline void printStats() {
		std::cout << "Performance Stats" << std::endl;
		for (auto& t : timers) {
			std::cout << t.first << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(t.second).count() << " ms, " << calls[t.first] << std::endl;
		}
	}
};

#define PROFILE 0
#if PROFILE
#define PROFILE_SCOPE(name) Recorder::Timer tt(name)
#else
#define PROFILE_SCOPE(name)
#endif
