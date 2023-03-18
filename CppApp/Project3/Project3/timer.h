#pragma once
#include <chrono>

class Timer
{
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_StartTime;

public:
    inline void Start();
    inline float GetDuration();
};


inline void Timer::Start()
{
    m_StartTime = std::chrono::high_resolution_clock::now();
}

inline float Timer::GetDuration()
{
    auto duration = std::chrono::high_resolution_clock::now() - m_StartTime;
    return duration.count();
}