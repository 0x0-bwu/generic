/**
 * @file Tools.hpp
 * @author bwu
 * @brief Some tools
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Macros.hpp"
#include "Units.hpp"

#if defined(GENERIC_OS_LINUX) || defined(GENERIC_OS_MAC)
#include <sys/stat.h>
#include <time.h>
#endif
#include <iostream>
#include <thread>
#include <chrono>
#include <atomic>
#include <ctime>
namespace generic {
///@brief tools and basic utilities
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

///@brief Local Time
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

///@brief Greenwich Mean Time
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

///@brief get current thread id
inline size_t ThreadId() noexcept
{
    return static_cast<size_t>(std::hash<std::thread::id>()(std::this_thread::get_id()));
}

///@brief sleep specified millisecond in current thread
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

///@brief timer counter from this object construct to destroy, usually used to track process execute time
class ProgressTimer
{
public:
    ///@brief constructs a prog\ess timer with specified out stream and display time unit
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

///@brief thread safe timer singleton, usually used to track total function execution time in threads
class AccumulatedTimer
{
public:
    class SubTimer
    {
    public:
        explicit SubTimer(std::shared_ptr<AccumulatedTimer> master) : m_master(master)
        {
            m_start = std::chrono::steady_clock::now(); 
        }
        
        ~SubTimer() { try { Stop(); } catch (...) {} }

        ///@brief stop timing explicitly, or it will stop when destructing
        void Stop()
        {
            if (not m_stop) {
                m_stop = true;
                auto end = std::chrono::steady_clock::now();
                auto elapse = std::chrono::duration_cast<std::chrono::nanoseconds>(end - m_start);
                m_master->AccumulateImp(elapse.count());
            }
        }

    private:
        bool m_stop = false;
        std::shared_ptr<AccumulatedTimer> m_master = nullptr;
        std::chrono::time_point<std::chrono::steady_clock> m_start;
    };

    friend SubTimer;
    AccumulatedTimer(const AccumulatedTimer &) = delete;
    AccumulatedTimer & operator= (const AccumulatedTimer &) = delete;
    
    ///@brief set time unit of returning timing count
    static void SetUnit(unit::Time unit)
    {
        Instance()->SetUnitImp(unit);
    }

    ///@brief get time unit of returning timing count
    static unit::Time GetUnit()
    {
        return Instance()->GetUnitImp();
    }

    ///@brief insert a timer object and start timing immediately
    ///@note should use a lvalue to hold the returning value explicitly, otherwise the timer will not work
    static std::unique_ptr<SubTimer> InsertTimer()
    {    
        return std::make_unique<SubTimer>(Instance());
    }

    ///@brief returning current timing count
    static double Count()
    {
        return Instance()->CountImp();
    }
    
    ///@brief reset timing count to zero
    static void Reset()
    {
        Instance()->ResetImp();
    }

private:
    AccumulatedTimer() = default;

    static std::shared_ptr<AccumulatedTimer> Instance()
    {
        static std::shared_ptr<AccumulatedTimer> timer(new AccumulatedTimer);
        return timer;
    }

    void SetUnitImp(unit::Time unit)
    {
        m_unit = unit;
    }

    unit::Time GetUnitImp() const
    {
        return m_unit;
    }

    void ResetImp()
    {
        m_count.store(0);
    }

    void AccumulateImp(int64_t time)
    {
        m_count.fetch_add(time);
    }

    double CountImp() const
    {
        double count = m_count.load() * unit::Scale2Second(unit::Time::Nanosecond);
        double scale = 1.0 / unit::Scale2Second(m_unit);
        return scale * count;
    }

private:
    unit::Time m_unit{unit::Time::Second};
    std::atomic<int64_t> m_count{0};//nanoseconds;
};

}//namespace tools
}//namespace generic