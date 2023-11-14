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

#ifndef BOOST_CHRONO_HEADER_ONLY
#define BOOST_CHRONO_HEADER_ONLY
#endif//BOOST_CHRONO_HEADER_ONLY
#include <boost/chrono/thread_clock.hpp>
#include <iostream>
#include <memory>
#include <thread>
#include <chrono>
#include <atomic>
#include <ctime>
namespace generic {
///@brief tools and basic utilities
namespace tools {

struct SystemClock
{
    using TimePoint = std::chrono::time_point<
                        std::chrono::system_clock, 
                            std::chrono::system_clock::duration>;
    static TimePoint Now() noexcept
    {
        return std::chrono::system_clock::now();
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
    ProgressTimer(std::string label = {}, unit::Time displayUnit = unit::Time::Second, std::ostream & os = std::cout)
     : m_label(label.empty() ? "progress" : label), m_os(os), m_timer(displayUnit)
    {
    }

    ~ProgressTimer()
    {
        try {
            m_os << m_label << " time: " << m_timer.Count();
            m_os << ::toString(m_timer.Unit());
            m_os << std::endl;
        }
        catch (...) {}
    }

private:
    std::string m_label;
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
            m_wtStart = std::chrono::steady_clock::now();
            m_ctStart = boost::chrono::thread_clock::now(); 
        }
        
        ~SubTimer() { try { Stop(); } catch (...) {} }

        ///@brief stop timing explicitly, or it will stop when destructing
        void Stop()
        {
            if (not m_stop) {
                m_stop = true;
                auto wtElapse = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - m_wtStart);
                auto ctElapse = boost::chrono::duration_cast<boost::chrono::nanoseconds>(boost::chrono::thread_clock::now() - m_ctStart);
                m_master->AccumulateImp(wtElapse.count(), ctElapse.count());
            }
        }

    private:
        bool m_stop = false;
        std::shared_ptr<AccumulatedTimer> m_master = nullptr;
        std::chrono::time_point<std::chrono::steady_clock> m_wtStart;
        boost::chrono::thread_clock::time_point m_ctStart;
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

    ///@brief returning current wall time count
    static double WallTime()
    {
        return Instance()->WallTimeImp();
    }

    ///@brief returning current cput time count
    static double CpuTime()
    {
        return Instance()->CpuTimeImp();
    }

    ///@brief returning current timing count [wall time, cpu time]
    static std::pair<double, double> Count()
    {
        return { WallTime(), CpuTime() };
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
        m_wtCount.store(0);
        m_ctCount.store(0);
    }

    void AccumulateImp(int64_t wallTime, int64_t cpuTime)
    {
        m_wtCount.fetch_add(wallTime);
        m_ctCount.fetch_add(cpuTime);
    }

    double WallTimeImp() const
    {
        double count = m_wtCount.load() * unit::Scale2Second(unit::Time::Nanosecond);
        double scale = 1.0 / unit::Scale2Second(m_unit);
        return scale * count;
    }

    double CpuTimeImp() const
    {
        double count = m_ctCount.load() * unit::Scale2Second(unit::Time::Nanosecond);
        double scale = 1.0 / unit::Scale2Second(m_unit);
        return scale * count;
    }

private:
    unit::Time m_unit{unit::Time::Second};
    std::atomic<int64_t> m_wtCount{0};//nanoseconds;
    std::atomic<int64_t> m_ctCount{0};//nanoseconds;
};

}//namespace tools
}//namespace generic