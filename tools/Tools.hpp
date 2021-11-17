#ifndef GENERIC_TOOLS_HPP
#define GENERIC_TOOLS_HPP
#include "common/Macros.hpp"
#include "Units.hpp"

#if defined(GENERIC_OS_LINUX) || defined(GENERIC_OS_MAC)
#include <sys/stat.h>
#include <time.h>
#endif

#include <iostream>
#include <thread>
#include <chrono>
#include <ctime>
namespace generic {
namespace tools {

struct SystemClock
{
    using Clock = typename std::chrono::system_clock;
    using Duration = typename Clock::duration;
    using Rep = typename Duration::rep;
    using Period = typename Duration::period;
    using TimePoint = std::chrono::time_point<Clock, Duration>;

    static TimePoint Now() noexcept
    {
        return Clock::now();
    }
};

inline std::tm LocalTime(const std::time_t & time)
{
    std::tm t;
#ifdef GENERIC_WINDOWS
    ::localtime_s(&t, &time);
#else
    ::localtime_r(&time, &t);
#endif
    return t;
}

//Greenwich Mean Time
inline std::tm GMT(const std::time_t & time)
{
    std::tm t;
#ifdef GENERIC_WINDOWS
    ::gmtime_s(&t, &time);
#else
    ::gmtime_r(&time, &t);
#endif
    return t;
}

inline size_t ThreadId() noexcept
{
    return static_cast<size_t>(std::hash<std::thread::id>()(std::this_thread::get_id()));
}

inline void SleepMilliseconds(size_t time)
{
    std::this_thread::sleep_for(std::chrono::milliseconds(time));
}

class Timer
{
public:
    explicit Timer(unit::Time unit = unit::Time::Second) : m_unit(unit) 
    {
        m_start = std::chrono::steady_clock::now();
    }
    
    double Count() const
    {
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapse = end - m_start;
        double scale = 1.0 / unit::Scale2Second(m_unit);
        return elapse.count() * scale;
    }
    void Restart() { m_start = std::chrono::steady_clock::now(); }
    unit::Time Unit() const { return m_unit; }

private:
    std::chrono::time_point<std::chrono::steady_clock> m_start;
    unit::Time m_unit;
};

class ProgressTimer
{
public:
    ProgressTimer(unit::Time displayUnit = unit::Time::Second, std::ostream & os = std::cout)
     : m_os(os), m_timer(displayUnit)
    {
    }

    ~ProgressTimer()
    {
        try {
            m_os << "progress time: " << m_timer.Count();
            m_os << ::toString(m_timer.Unit());
            m_os << std::endl;
        }
        catch (...) {}
    }

private:
    std::ostream & m_os;
    Timer m_timer;
};

}//namespace tools
}//namespace generic
#endif//GENERIC_TOOLS_HPP