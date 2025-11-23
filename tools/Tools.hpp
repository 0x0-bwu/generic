/**
 * @file Tools.hpp
 * @author bwu
 * @brief Some tools
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Exception.hpp"
#include "Units.hpp"

#if defined(GENERIC_OS_LINUX) || defined(GENERIC_OS_MAC)
#include <sys/stat.h>
#include <time.h>
#endif

#include <iostream>
#include <vector>
#include <thread>
#include <memory>
#include <chrono>
#include <atomic>
#include <ctime>
namespace generic::tools {
///@brief tools and basic utilities
struct SystemClock
{
    using Clock = std::chrono::system_clock;
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
    explicit Timer(unit::Time unit = unit::Time::SECOND) : m_unit(unit) 
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
/**
 * @brief Brief description of Restart.
 * @return void
 */
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
    ProgressTimer(std::string label = {}, unit::Time displayUnit = unit::Time::SECOND, std::ostream & os = std::cout)
     : m_label(label.empty() ? "progress" : label), m_os(os), m_timer(displayUnit)
    {
    }

    ~ProgressTimer()
    {
        try {
            m_os << m_label << " time: " << m_timer.Count();
            m_os << generic::toString(m_timer.Unit());
            m_os << std::endl;
        }
        catch (...) {}
    }

private:
    std::string m_label;
    std::ostream & m_os;
    Timer m_timer;
};

struct AccumulatedTimerTag
{
    inline static constexpr size_t groups = 1;
};
///@brief thread safe timer singleton, usually used to track total function execution time in threads
template <typename Tag = AccumulatedTimerTag>
class AccumulatedTimer
{
public:
    struct ExampleTag
    {
        inline static constexpr size_t groups = 1;
    };

    static_assert(Tag::groups > 0, "You should define a constexpr size_t type of \"groups\" which at least equals 1 in your tag type");
    class SubTimer
    {
    public:
        explicit SubTimer(std::shared_ptr<AccumulatedTimer<Tag>> master, size_t group)
         : m_group(group), m_master(master)
        {
            m_start = std::chrono::steady_clock::now();
        }
        
        ~SubTimer() { try { Stop(); } catch (...) {} }

        ///@brief stop timing explicitly, or it will stop when destructing
        void Stop()
        {
            if (not m_stop) {
                m_stop = true;
                auto elapse = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - m_start);
                m_master->AccumulateImp(m_group, elapse.count());
            }
        }

    private:
        bool m_stop = false;
        std::size_t m_group = 0;
        std::chrono::time_point<std::chrono::steady_clock> m_start;
        std::shared_ptr<AccumulatedTimer<Tag>> m_master = nullptr;
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
    ///@param group, specify the timer group, result will be added to corresponding timer group
    ///@note should use a lvalue to hold the returning value explicitly, otherwise the timer will not work
    static std::unique_ptr<SubTimer> InsertTimer(size_t group = 0)
    {    
        return std::make_unique<SubTimer>(Instance(), group);
    }

    ///@brief return current timing count [wall times, call times] of specified group
    static std::pair<double, size_t> Count(size_t group = 0)
    {
        return Instance()->CountImp(group);
    }
    
    ///@brief reset specified group timing count to zero
    static void Reset(size_t group)
    {
        Instance()->ResetImp(group);
    }

    ///@brief reset all groups timing count to zero
    static void Reset()
    {
        Instance()->ResetImp();
    }

private:
    AccumulatedTimer() { ResetImp(); }

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

    void ResetImp(size_t group)
    {
        m_times[group].store(0);
        m_count[group].store(0);
    }

    void ResetImp()
    {
        for (auto & t : m_times) t.store(0);
        for (auto & c : m_count) c.store(0);
    }

    void AccumulateImp(size_t group, int64_t time)
    {
        m_times[group] += time;
        m_count[group] += 1;
    }

    std::pair<double, size_t> CountImp(size_t group) const
    {
        double times = m_times.at(group) * unit::Scale2Second(unit::Time::NANOSECOND);
        double scale = 1.0 / unit::Scale2Second(m_unit);
        return { scale * times, m_count.at(group)};
    }

private:
    unit::Time m_unit{unit::Time::SECOND};
    std::array<std::atomic<int64_t>, Tag::groups> m_count;
    std::array<std::atomic<int64_t>, Tag::groups> m_times;//nanoseconds;
};

using ThreadTimer = AccumulatedTimer<AccumulatedTimerTag>;

}//namespace generic::tools
