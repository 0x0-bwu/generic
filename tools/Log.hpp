/**
 * @file Log.hpp
 * @author bwu
 * @brief A header only log library
 * @version 0.1
 * @date 2022-02-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include "generic/common/Exception.hpp"
#include "FileSystem.hpp"
#include "Format.hpp"
#include "Tools.hpp"
#include <unordered_map>
#include <functional>
#include <cstring>
#include <cctype>
#include <memory>
#include <vector>
#include <string>
#include <atomic>
#include <cstdio>
#include <ctime>
#include <mutex>

#define GENERIC_LOG_LEVEL_TRACE 0
#define GENERIC_LOG_LEVEL_DEBUG 1
#define GENERIC_LOG_LEVEL_INFO  2
#define GENERIC_LOG_LEVEL_WARN  3
#define GENERIC_LOG_LEVEL_ERROR 4
#define GENERIC_LOG_LEVEL_FATAL 5
#define GENERIC_LOG_LEVEL_OFF   6

namespace generic {
///@brief namespace of log library

namespace log {

enum class Level
{
    Trace = GENERIC_LOG_LEVEL_TRACE,
    Debug = GENERIC_LOG_LEVEL_DEBUG,
    Info  = GENERIC_LOG_LEVEL_INFO,
    Warn  = GENERIC_LOG_LEVEL_WARN,
    Error = GENERIC_LOG_LEVEL_ERROR,
    Fatal = GENERIC_LOG_LEVEL_FATAL,
    Off   = GENERIC_LOG_LEVEL_OFF,
    nLevels
};

} // namespace log

inline std::string toString(log::Level level)
{
    using namespace log;
    switch (level)
    {
        case Level::Trace : return "Trace";
        case Level::Debug : return "Debug";
        case Level::Info  : return "Info" ;
        case Level::Warn  : return "Warn" ;
        case Level::Error : return "Error";
        case Level::Fatal : return "Fatal";
        case Level::Off   : return "Off"  ;
        default: return "";
    }
}
namespace log {

using LogClock = tools::SystemClock;
using LogTimePoint = typename LogClock::TimePoint;
using LogLevel = std::atomic<int>;
using ErrHandler = std::function<void(const std::string & err)>;
using Format = boost::format;

namespace details {

struct SourceLoc
{
    SourceLoc() = default;
    SourceLoc(const char * _file, int _line, const char * _func)
     : line(_line), file(_file), func(_func) {}
    
    bool Empty() const noexcept { return 0 == line; }

    int line{0};
    const char * file{nullptr};
    const char * func{nullptr};
};

struct LogMsg
{    
    LogMsg() = default;

    LogMsg(std::string logger, Level level, std::string msg)
     : LogMsg(SourceLoc{}, logger, level, msg)
    {}

    LogMsg(SourceLoc source, std::string logger, Level level, std::string msg)
     : LogMsg(LogClock::Now(), source, logger, level, msg)
    {}

    LogMsg(LogTimePoint _time, SourceLoc _source, std::string _logger, Level _level, std::string msg)
     : logger(_logger), level(_level), time(_time), threadId(tools::ThreadId()), source(_source), payload(msg)
    {}
    
    std::string logger;
    Level level = Level::Off;
    LogTimePoint time;
    size_t threadId = 0;
    SourceLoc source;
    std::string payload;
};

class Formatter
{
public:
    virtual ~Formatter() = default;
    virtual void Format(const LogMsg & msg, std::string & out) = 0;
    virtual std::unique_ptr<Formatter> Clone() const = 0;
};

struct PaddingInfo
{
    enum class PadSide { Left, Right, Center };
    PaddingInfo() = default;
    PaddingInfo(size_t _width, PadSide _side, bool _truncate)
     : width(_width), side(_side), truncate(_truncate), enabled(true)
    {}

    size_t width = 0;
    PadSide side = PadSide::Left;
    bool truncate = false;
    bool enabled = false;
};

class ScopedPadder
{
public:
    ScopedPadder(size_t wrappedSize, const PaddingInfo & padInfo, std::string & dest)
     : m_padInfo(padInfo), m_dest(dest)
    {
        m_remainingPad = static_cast<long>(padInfo.width) - static_cast<long>(wrappedSize);
        if(m_remainingPad <= 0) return;

        if(m_padInfo.side == PaddingInfo::PadSide::Left){
            PadIt(m_remainingPad);
            m_remainingPad = 0;
        }
/**
 * @brief Brief description of if.
 * @param PaddingInfo::PadSide::Center
 * @return else
 */
        else if(m_padInfo.side == PaddingInfo::PadSide::Center){
            auto halfPad = m_remainingPad / 2;
            auto reminder = m_remainingPad & 1;
            PadIt(halfPad);
            m_remainingPad = halfPad + reminder;
        }
    }

private:
    void PadIt(long count)
    {
        m_dest.append(m_spaces.substr(0, count));
    }
private:
    const PaddingInfo & m_padInfo;
    std::string & m_dest;
    long m_remainingPad;
    std::string m_spaces = std::string(64, ' ');
};

struct NullScopedPadder
{
   NullScopedPadder(size_t, const PaddingInfo & , std::string & ) {}

   template <typename T>
/**
 * @brief Brief description of CountDigits.
 * @param T
 * @return static unsigned int
 */
   static unsigned int CountDigits(T) { return 0; }
};

class FlagFormatter
{
public:
    FlagFormatter() = default;
/**
 * @brief Brief description of FlagFormatter.
 * @param m_paddingInfo(info
 * @return explicit
 */
    explicit FlagFormatter(PaddingInfo info) : m_paddingInfo(info) {}
    virtual ~FlagFormatter() = default;
    virtual void Format(const LogMsg & msg, const std::tm & time, std::string & dest) = 0;
protected:
    PaddingInfo m_paddingInfo;
};

class AggregateFormatter final : public FlagFormatter
{
public:
    AggregateFormatter() = default;

/**
 * @brief Brief description of AddChar.
 * @param c
 * @return void
 */
    void AddChar(char c) { m_str.push_back(c); }

    void Format(const LogMsg & msg, const std::tm & , std::string & dest) override
    {
        dest.append(m_str);
    }
private:
    std::string m_str;
};

class FullFormatter final : public FlagFormatter
{
    using FormatHelper = fmt::FormatHelper;
public:
/**
 * @brief Brief description of FullFormatter.
 * @param FlagFormatter(info
 * @return explicit
 */
    explicit FullFormatter(PaddingInfo info) : FlagFormatter(info) {}

    void Format(const LogMsg & msg, const std::tm & time, std::string & dest) override
    {
        using std::chrono::duration_cast;
        using std::chrono::microseconds;
        using std::chrono::seconds;

        auto millisFraction = [](LogTimePoint tp)
        {
            auto duration = tp.time_since_epoch();
            auto secs = duration_cast<seconds>(duration);
            return duration_cast<microseconds>(duration) - duration_cast<microseconds>(secs);
        };

        auto duration = msg.time.time_since_epoch();
        auto secs = duration_cast<seconds>(duration);
        
        if(m_cacheTimeStamp != secs || m_cacheDataTime.size() == 0){
            m_cacheDataTime.clear();
            m_cacheDataTime.push_back('[');
            FormatHelper::AppendInt<4>(m_cacheDataTime, time.tm_year + 1900);
            m_cacheDataTime.push_back('-');
            FormatHelper::AppendInt<2>(m_cacheDataTime, time.tm_mon + 1, '0');
            m_cacheDataTime.push_back('-');
            FormatHelper::AppendInt<2>(m_cacheDataTime, time.tm_mday, '0');
            m_cacheDataTime.push_back(' ');
            FormatHelper::AppendInt<2>(m_cacheDataTime, time.tm_hour, '0');
            m_cacheDataTime.push_back(':');
            FormatHelper::AppendInt<2>(m_cacheDataTime, time.tm_min, '0');
            m_cacheDataTime.push_back(':');
            FormatHelper::AppendInt<2>(m_cacheDataTime, time.tm_sec, '0');
            m_cacheDataTime.push_back('.');

            m_cacheTimeStamp = secs;
        }
        dest.append(m_cacheDataTime.begin(), m_cacheDataTime.end());
        
        auto millis = millisFraction(msg.time);
        FormatHelper::AppendInt<3>(dest, millis.count(), '0');
        dest.push_back(']');
        dest.push_back(' ');

        if(!msg.logger.empty()){
            dest.push_back('[');
            dest.append(msg.logger);
            dest.push_back(']');
            dest.push_back(' ');
        }

        dest.push_back('[');
        dest.append(toString(msg.level));
        dest.push_back(']');
        dest.push_back(' ');

        if(!msg.source.Empty()){
            dest.push_back('[');
            dest.append(fs::FileName(msg.source.file));
            dest.push_back(']');
            dest.push_back(' ');
        }

        dest.append(msg.payload);
    }

private:
    std::chrono::seconds m_cacheTimeStamp{0};
    std::string m_cacheDataTime;
};

class CustomFlagFormatter : public FlagFormatter
{
public:
    virtual std::unique_ptr<CustomFlagFormatter> Clone() const = 0;
/**
 * @brief Brief description of SetPaddingInfo.
 * @param info
 * @return void
 */
    void SetPaddingInfo(PaddingInfo info) { FlagFormatter::m_paddingInfo = info; }
};

enum class PatternTimeType { Local, UTC };

class PatternFormatter : public Formatter
{
public:
    using CustomFlags = std::unordered_map<char, std::unique_ptr<CustomFlagFormatter> >;
    explicit PatternFormatter(std::string pattern, PatternTimeType timeType = PatternTimeType::Local,
                              std::string eol = GENERIC_DEFAULT_EOL, CustomFlags customFlags = CustomFlags{})
     : m_pattern(pattern), m_eol(std::move(eol)), m_timeType(timeType), m_lastLogSec(0), m_customHandlers(std::move(customFlags))
    {
        std::memset(&m_cachedTime, 0, sizeof(m_cachedTime));
        CompilePattern(m_pattern);
    }

    explicit PatternFormatter(PatternTimeType timeType = PatternTimeType::Local, std::string eol = GENERIC_DEFAULT_EOL)
     : m_pattern("%+"), m_eol(std::move(eol)), m_timeType(timeType), m_lastLogSec(0)
    {
         std::memset(&m_cachedTime, 0, sizeof(m_cachedTime));
         m_formatters.push_back(std::unique_ptr<FlagFormatter>(new FullFormatter(PaddingInfo{})));
    }

    std::unique_ptr<Formatter> Clone() const
    {
        CustomFlags clonedCumstomFlags;
        for(auto & it : m_customHandlers){
            clonedCumstomFlags[it.first] = it.second->Clone();
        }
        return std::unique_ptr<Formatter>(new PatternFormatter(m_pattern, m_timeType, m_eol, std::move(clonedCumstomFlags)));
    }

    void Format(const LogMsg & msg, std::string & dest)
    {
        auto secs = std::chrono::duration_cast<std::chrono::seconds>(msg.time.time_since_epoch());
        if(secs != m_lastLogSec){
            m_cachedTime = GetTime(msg);
            m_lastLogSec = secs;
        }

        for(auto & f : m_formatters){
            f->Format(msg, m_cachedTime, dest);
        }
        dest.append(GENERIC_DEFAULT_EOL);
    }
    
    void CompilePattern(const std::string & pattern)
    {
        auto end = pattern.end();
        std::unique_ptr<AggregateFormatter> userChars;
        m_formatters.clear();
        for(auto iter = pattern.begin(); iter != end; ++iter){
            if(*iter == '%'){
                if(userChars){
                    m_formatters.push_back(std::move(userChars));
                }
                
                auto padding = HandlePadSpec(++iter, end);

                if(iter != end){
                    if(padding.enabled){
                        HandleFlag<ScopedPadder>(*iter, padding);
                    }
                    else HandleFlag<NullScopedPadder>(*iter, padding);
                }
                else break;
            }
            else {
                if(!userChars)
                    userChars.reset(new AggregateFormatter);
                userChars->AddChar(*iter);
            }
        }
        if(userChars)
            m_formatters.push_back(std::move(userChars));
    }

private:
    std::tm GetTime(const LogMsg & msg)
    {
        if(m_timeType == PatternTimeType::Local){
            return tools::LocalTime(LogClock::Clock::to_time_t(msg.time));
        }
        return tools::GMT(LogClock::Clock::to_time_t(msg.time));
    }

    PaddingInfo HandlePadSpec(std::string::const_iterator & iter, std::string::const_iterator end)
    {
        const size_t maxWidth = 64;
        if(iter == end) return PaddingInfo{};

        PaddingInfo::PadSide side;
        switch(*iter){
            case '-' :
                side = PaddingInfo::PadSide::Right;
                ++iter;
                break;
            case '=' :
                side = PaddingInfo::PadSide::Center;
                ++iter;
                break;
            default:
                side = PaddingInfo::PadSide::Left;
                break;
        }

        if(iter == end || !std::isdigit(static_cast<unsigned char>(*iter)))
            return PaddingInfo{};
        
        auto width = static_cast<size_t>(*iter) - '0';
        for(++iter; iter != end && std::isdigit(static_cast<unsigned char>(*iter)); ++iter){
            auto digit = static_cast<size_t>(*iter) - '0';
            width = width * 10 + digit;
        }

        bool truncate;
        if(iter != end && *iter == '!'){
            truncate = true;
            ++iter;
        }
        else truncate = false;

        return PaddingInfo {std::min<size_t>(width, maxWidth), side, truncate};
    }

    template <typename Padder>
    void HandleFlag(char flag, PaddingInfo padding)
    {
       auto iter = m_customHandlers.find(flag);
       if(iter != m_customHandlers.end()){
           auto customHandler = iter->second->Clone();
           customHandler->SetPaddingInfo(padding);
           m_formatters.push_back(std::move(customHandler));
       }

       switch (flag)
       {
       case ('+') :
            m_formatters.push_back(std::unique_ptr<FlagFormatter>(new FullFormatter(padding)));
            break;
            //todo  
       default:
            break;
       }
    }
private:
    std::string m_pattern;
    std::string m_eol;
    PatternTimeType m_timeType;
    std::tm m_cachedTime;
    std::chrono::seconds m_lastLogSec;
    std::vector<std::unique_ptr<FlagFormatter> > m_formatters;
    CustomFlags m_customHandlers;
};

template <typename T>
class CircularQueue
{
    size_t m_maxItems;
    size_t m_head = 0;
    size_t m_tail = 0;
    size_t m_overrunCounter = 0;
    std::vector<T> m_values;
public:
    using value_type = T;

    CircularQueue() = default;
/**
 * @brief Brief description of CircularQueue.
 * @param 1)
 * @param m_values(maxItems
 * @return explicit
 */
    explicit CircularQueue(size_t maxItems) : m_maxItems(maxItems + 1), m_values(maxItems) {}
    CircularQueue(const CircularQueue & other) = default;
    CircularQueue & operator= (const CircularQueue & other) = default;

    CircularQueue(CircularQueue && other) { Move(std::forward<CircularQueue>(other)); }
    CircularQueue & operator= (CircularQueue && other) { Move(std::forward<CircularQueue>(other)); return *this; }

    void PushBack(T && item)
    {
        if(m_maxItems > 0){
            m_values[m_tail] = std::move(item);
            m_tail = (m_tail + 1) & m_maxItems;
            if(m_tail == m_head){
                m_head = (m_head + 1) % m_maxItems;
                ++m_overrunCounter;
            }
        }
    }

/**
 * @brief Brief description of Front.
 * @return T &
 */
    T & Front() { return m_values[m_head]; }
    const T & Front() const { return m_values[m_head]; }
/**
 * @brief Brief description of Size.
 * @return size_t
 */
    size_t Size() { return m_tail > m_head ? (m_tail - m_head) : m_maxItems - (m_head - m_tail); }
    const T & At(size_t i) const
    {
        GENERIC_ASSERT(i < Size());
        return m_values[(m_head + i) % m_maxItems];
    }

/**
 * @brief Brief description of PopFront.
 * @return void
 */
    void PopFront() { m_head = (m_head + 1) % m_maxItems; }
    bool Empty() const { return m_tail == m_head; }
    bool Full() const
    {
        if(m_maxItems > 0) return ((m_tail + 1) % m_maxItems) == m_head;
        return false;
    }

    size_t OverrunCounter() const { return m_overrunCounter; }

private:
    void Move(CircularQueue && other)
    {
        m_maxItems = other.m_maxItems;
        m_head = other.m_head;
        m_tail = other.m_tail;
        m_overrunCounter = other.m_overrunCounter;
        m_values = std::move(other.m_values);

        other.m_maxItems = 0;
        other.m_head = other.m_tail = 0;
        other.m_overrunCounter = 0;
    }
};

class BackTracer
{
    mutable std::mutex m_mutex;
    std::atomic_bool m_enabled;
    CircularQueue<LogMsg> m_messages;
public:
    BackTracer() : m_enabled(false) {}
    
    BackTracer(const BackTracer & other)
    {
        std::lock_guard<std::mutex> lock(other.m_mutex);
        m_enabled = other.Enabled();
        m_messages = other.m_messages;
    }

    BackTracer(BackTracer && other)
    {
        Move(std::forward<BackTracer>(other));
    }

    BackTracer & operator=(BackTracer other)
    {
        Move(std::move(other));
        return *this;
    }

    void Enable(size_t size)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_enabled.store(true, std::memory_order_relaxed);
        m_messages = CircularQueue<LogMsg>{size};
    }

    void Disable()
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_enabled.store(false, std::memory_order_relaxed);
    }

    bool Enabled() const
    { 
        return m_enabled.load(std::memory_order_relaxed);
    }

    void PushBack(LogMsg && msg)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_messages.PushBack(std::forward<LogMsg>(msg));
    }

    void PushBack(const LogMsg & msg)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_messages.PushBack(LogMsg(msg));
    }

    void ForeachPop(std::function<void(const LogMsg &)> func)
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        while(!m_messages.Empty()){
            func(m_messages.Front());
            m_messages.PopFront();
        }
    }

private:
    void Move(BackTracer && other)
    {
        std::lock_guard<std::mutex> lock(other.m_mutex);
        m_enabled = other.Enabled();
        m_messages = std::move(other.m_messages);
    }
};

}//namespace details

///@brief skink definitions
namespace sinks {

class Sink
{
public:
    virtual ~Sink() = default;
    virtual void Log(const details::LogMsg & msg) = 0;
    virtual void Flush() = 0;
    virtual void SetPattern(const std::string & pattern) = 0;
    virtual void SetFormatter(std::unique_ptr<details::Formatter> formatter) = 0;

    void SetLevel(Level level)
    {
        m_level.store(static_cast<int>(level), std::memory_order_relaxed);
    }

    Level GetLevel() const
    {
        return static_cast<Level>(m_level.load(std::memory_order_relaxed));
    }

    bool ShouldLog(Level level) const
    {
        return static_cast<int>(level) >= m_level.load(std::memory_order_relaxed);
    }

protected:
    LogLevel m_level{static_cast<int>(Level::Trace)};
};

template <typename mutex_t>
class BaseSink : public Sink
{
    using LogMsg = details::LogMsg;
public:
    BaseSink() : m_formatter(new details::PatternFormatter) {}
/**
 * @brief Brief description of BaseSink.
 * @param m_formatter(std::move(formatter)
 * @return explicit
 */
    explicit BaseSink(std::unique_ptr<details::Formatter> formatter) : m_formatter(std::move(formatter)) {}
    ~BaseSink() override = default;

    BaseSink(const BaseSink & ) = delete;
    BaseSink(BaseSink && ) = delete;

    BaseSink & operator=(const BaseSink & ) = delete;
    BaseSink & operator=(BaseSink && ) = delete;

    void Log(const LogMsg & msg) final
    {
        std::lock_guard<mutex_t> lock(m_mutex);
        SinkImp(msg);
    }

    void Flush() final
    {
        std::lock_guard<mutex_t> lock(m_mutex);
        FlushImp();
    }

    void SetPattern(const std::string & pattern) final
    {
        std::lock_guard<mutex_t> lock(m_mutex);
        SetPatternImp(pattern);
    }

    void SetFormatter(std::unique_ptr<details::Formatter> formatter) final
    {
        std::lock_guard<mutex_t> lock(m_mutex);
        SetFormatterImp(std::move(formatter));
    }


protected:
    virtual void SinkImp(const LogMsg & msg) = 0;
    virtual void FlushImp() = 0;
    virtual void SetPatternImp(const std::string & pattern)
    {
        SetFormatterImp(std::unique_ptr<details::Formatter>(new details::PatternFormatter(pattern)));
    }
    virtual void SetFormatterImp(std::unique_ptr<details::Formatter> fomatter)
    {
        m_formatter = std::move(fomatter);
    }

protected:
    std::unique_ptr<details::Formatter> m_formatter;
    mutable mutex_t m_mutex;
};

template <typename mutex_t>
class BasicFileSink final : public BaseSink<mutex_t>
{
    using LogMsg = details::LogMsg;
public:
    explicit BasicFileSink(std::string_view filename, bool truncate = true)
    {
        m_fileHelper.Open(filename, truncate);
    }

    const std::string & Filename() const { return m_fileHelper.Filename(); }

protected:
    void SinkImp(const LogMsg & msg) override
    {
        std::string formatted;
        BaseSink<mutex_t>::m_formatter->Format(msg, formatted);
        m_fileHelper.Write(formatted);
    }
    void FlushImp() override
    {
        m_fileHelper.Flush();
    }
private:
    fs::FileHelper m_fileHelper;
};

template <typename mutex_t>
class BasicOstreamSink final : public BaseSink<mutex_t>
{
    using LogMsg = details::LogMsg;
public:
    explicit BasicOstreamSink(std::ostream & os, bool forceFlush = false)
     : m_ostream(os), m_forceFlush(forceFlush)
    {}

    BasicOstreamSink(const BasicOstreamSink & ) = delete;
    BasicOstreamSink & operator= (const BasicOstreamSink & ) = delete;

protected:
    void SinkImp(const LogMsg & msg) override
    {
        std::string formatted;
        BaseSink<mutex_t>::m_formatter->Format(msg, formatted);
        m_ostream.write(formatted.data(), static_cast<std::streamsize>(formatted.size()));
        if(m_forceFlush) m_ostream.flush();
    }

    void FlushImp() override
    {
        m_ostream.flush();
    }

    std::ostream & m_ostream;
    bool m_forceFlush;
};

}//namespace sinks

using SinkPtr = std::shared_ptr<sinks::Sink>;
using SinksInitializerList = std::initializer_list<SinkPtr>;

class Logger
{
    using LogMsg = details::LogMsg;
    using SourceLoc = details::SourceLoc;
protected:
    std::string m_name;
    std::vector<SinkPtr> m_sinks;
    LogLevel m_level = static_cast<int>(Level::Trace);
    LogLevel m_flushLevel = static_cast<int>(Level::Off);
    details::BackTracer m_tracer;
    ErrHandler m_errHandler = nullptr;
public:
    ///@brief constructs an empty logger
    explicit Logger(std::string name) : m_name(std::move(name)), m_sinks() {}
    template <typename iterator>
    ///@brief constructs a logger with range on sinks
    Logger(std::string name, iterator begin, iterator end) : m_name(std::move(name)), m_sinks(begin, end) {}
    ///@brief constructs a logger with single sink
    Logger(std::string name, SinkPtr sink) : m_name(std::move(name)), m_sinks({std::move(sink)}) {}
    ///@brief constructs a logger with a list of sinks
    Logger(std::string name, SinksInitializerList sinks) : Logger(std::forward<std::string>(name), sinks.begin(), sinks.end()) {}
    virtual ~Logger() = default;

    void Trace(const std::string & msg)
    {
        return Log(Level::Trace, msg);
    }

    void Debug(const std::string & msg)
    {
        return Log(Level::Debug, msg);
    }

    void Info(const std::string & msg)
    {
        return Log(Level::Info, msg);
    }

    void Warn(const std::string & msg)
    {
        return Log(Level::Warn, msg);
    }

    void Error(const std::string & msg)
    {
        return Log(Level::Error, msg);
    }

    void Fatal(const std::string & msg)
    {
        return Log(Level::Fatal, msg);
    }
    
    void Log(Level level, const std::string & msg) 
    {
        return Log(SourceLoc{}, level, msg);
    }

    void Log(SourceLoc loc, Level level, std::string msg)
    {
        bool logEnabled = ShouldLog(level);
        bool tracebackEnabled = m_tracer.Enabled();
        if(!logEnabled && !tracebackEnabled) return;

        details::LogMsg logMsg(loc, m_name, level, msg);
        LogOne(logMsg, logEnabled, tracebackEnabled);
    }

    template <typename... Args>
    void Trace(const std::string & format, Args &&... args)
    {
        Log(Level::Trace, format, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Debug(const std::string & format, Args &&... args)
    {
        Log(Level::Debug, format, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Info(const std::string & format, Args &&... args)
    {
        Log(Level::Info, format, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Warn(const std::string & format, Args &&... args)
    {
        Log(Level::Warn, format, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Error(const std::string & format, Args &&... args)
    {
        Log(Level::Error, format, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Fatal(const std::string & format, Args &&... args)
    {
        Log(Level::Fatal, format, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Log(Level level, const std::string & format, Args &&... args)
    {
        Log(SourceLoc{}, level, format, std::forward<Args>(args)...);
    }

    template <typename... Args>
    void Log(SourceLoc loc, Level level, const std::string & format, Args &&... args)
    {
        LogOne(loc, level, format, std::forward<Args>(args)...);
    }

    bool ShouldLog(Level level) const
    {
        return static_cast<int>(level) >= m_level.load(std::memory_order_relaxed);
    }

    bool ShouldBackTrace() const
    { 
        return m_tracer.Enabled();
    }

    void SetLevel(Level level)
    {
        m_level.store(static_cast<int>(level));
    }

    void FlushOn(Level level)
    {
        m_flushLevel.store(static_cast<int>(level));
    }

    void EnableBacktrace(size_t nMessages)
    {
        m_tracer.Enable(nMessages);
    }

    Level GetLevel() const
    {
        return static_cast<Level>(m_level.load(std::memory_order_relaxed));
    }

    std::string Name() const
    {
        return m_name;
    }

    void SetFormatter(std::unique_ptr<details::Formatter> formatter)
    {
        for(auto iter = m_sinks.begin(); iter != m_sinks.end(); ++iter ){
            if(std::next(iter) == m_sinks.end()){
                (*iter)->SetFormatter(std::move(formatter));
                break;
            }
            else (*iter)->SetFormatter(formatter->Clone());
        }
    }

    void SetErrorHandler(ErrHandler handler)
    {
        m_errHandler = std::move(handler);
    }

    void DumpBacktrace()
    {
        DumpBacktraceImp();
    }

protected:
    void LogOne(const LogMsg & msg, bool logEnabled, bool tracebackEnabled)
    {
        if(logEnabled) SinkOne(msg);
        if(tracebackEnabled) TraceOne(msg);
    }

    template <typename... Args>
    void LogOne(SourceLoc loc, Level level, const std::string & format, Args &&... args)
    {
        bool logEnabled = ShouldLog(level);
        bool tracebackEnabled = ShouldBackTrace();
        if(!logEnabled && !tracebackEnabled) return;

        GENERIC_TRY {
            auto str = fmt::Fmt2Str(format, std::forward<Args>(args)...);
            auto msg = details::LogMsg(loc, m_name, level, str);
            LogOne(msg, logEnabled, tracebackEnabled);
        }
        GENERIC_CATCH
    }

    void SinkOne(const LogMsg & msg)
    {
        for(auto & sink : m_sinks){
            if(sink->ShouldLog(msg.level)){
                GENERIC_TRY{
                    sink->Log(msg);
                }
                GENERIC_CATCH
            }
        }

        if(ShouldFlush(msg)) FlushSinks();
    }

    void TraceOne(const LogMsg & msg)
    {
        m_tracer.PushBack(msg);
    }

    void FlushSinks()
    {
        for(auto & sink : m_sinks){
            GENERIC_TRY {
                sink->Flush();
            }
            GENERIC_CATCH
        }
    }

    void BackendFlush()
    {
        for(auto & sink : m_sinks){
            GENERIC_TRY {
                sink->Flush();
            }
            GENERIC_CATCH
        }
    }

    bool ShouldFlush(const LogMsg & msg)
    {
        auto flushLevel = m_flushLevel.load(std::memory_order_relaxed);
        return (static_cast<int>(msg.level) >= flushLevel) && (msg.level != Level::Off);
    }

    void DumpBacktraceImp()
    {
        if(m_tracer.Enabled()){
            SinkOne(LogMsg{ Name(), Level::Info, "****************** Backtrace Start ******************"});
            m_tracer.ForeachPop([this](const LogMsg & msg){ this->SinkOne(msg); });
            SinkOne(LogMsg{ Name(), Level::Info, "****************** Backtrace End ********************"});
        }
    }
};

class Registry
{
    using LogLevels = std::unordered_map<std::string, Level>;
    Registry() : m_formatter(new details::PatternFormatter)
    {
        std::string name = "";
        m_defaultLogger = std::make_shared<Logger>(name);
        m_loggers[name] = m_defaultLogger;
    }

    ~Registry() = default;

public:
    Registry(const Registry & ) = delete;
    Registry & operator=(const Registry & ) = delete;

    static Registry & Instance()
    {
        static Registry registry;
        return registry;
    }

    bool RegisterLogger(std::shared_ptr<Logger> newLogger)
    {
        std::lock_guard<std::mutex> lock(m_loggerMapMutex);
        return RegisterOneLogger(std::move(newLogger));

    }

    void InitializeLogger(std::shared_ptr<Logger> newLogger)
    {
        std::lock_guard<std::mutex> lock(m_loggerMapMutex);

        newLogger->SetFormatter(m_formatter->Clone());

        if(m_errHandler){
            newLogger->SetErrorHandler(m_errHandler);
        }

        auto iter = m_logLevels.find(newLogger->Name());
        auto newLevel = iter != m_logLevels.end() ? iter->second : m_globalLevel;
        newLogger->SetLevel(newLevel);
        newLogger->FlushOn(m_flushLevel);

        if(m_backtraceN > 0){
            newLogger->EnableBacktrace(m_backtraceN);
        }

        if(m_autoRegistration){
            RegisterOneLogger(std::move(newLogger));
        }
    }

    std::shared_ptr<Logger> Get(const std::string & name)
    {
        std::lock_guard<std::mutex> lock(m_loggerMapMutex);
        auto found = m_loggers.find(name);
        return found == m_loggers.end() ? nullptr : found->second;
    }

    std::shared_ptr<Logger> DefaultLogger()
    {
        std::lock_guard<std::mutex> lock(m_loggerMapMutex);
        return m_defaultLogger;
    }

    ///@note cannot be used concurrently with SetDefaultLogger(...)
    Logger * GetDefaultRaw() { return m_defaultLogger ? m_defaultLogger.get() : nullptr; }

    void SetDefaultLogger(std::shared_ptr<Logger> newDefaultLogger)
    {
        std::lock_guard<std::mutex> lock(m_loggerMapMutex);
        if(nullptr != m_defaultLogger){
            m_loggers.erase(m_defaultLogger->Name());
        }
        if(nullptr != newDefaultLogger){
            m_loggers[newDefaultLogger->Name()] = newDefaultLogger;
        }
        m_defaultLogger = std::move(newDefaultLogger);
    }

    void SetFormatter(std::unique_ptr<details::Formatter> formatter)
    {
        if(nullptr == formatter) return;

        std::lock_guard<std::mutex> lock(m_loggerMapMutex);
        m_formatter = std::move(formatter);
        for(auto & logger : m_loggers){
            logger.second->SetFormatter(m_formatter->Clone());
        }
    }

    void Drop(const std::string & name)
    {
        std::lock_guard<std::mutex> lock(m_loggerMapMutex);
        m_loggers.erase(name);
        if(m_defaultLogger && m_defaultLogger->Name() == name){
            m_defaultLogger.reset();
        }
    }

    void DropAll()
    {
        std::lock_guard<std::mutex> lock(m_loggerMapMutex);
        m_loggers.clear();
        m_defaultLogger.reset();
    }

    void ShutDown()
    {
        //todo periodic flush
        DropAll();
        //todo thread pool
    }

private:
    bool RegisterOneLogger(std::shared_ptr<Logger> newLogger)
    {
        auto name = newLogger->Name();
        if(m_loggers.count(name)) return false;
        m_loggers[name] = std::move(newLogger);
        return true;
    }

private:
    mutable std::mutex m_loggerMapMutex;
    mutable std::mutex m_flusherMutex;
    std::recursive_mutex m_tmpMutex;
    std::unordered_map<std::string, std::shared_ptr<Logger> > m_loggers;
    std::shared_ptr<Logger> m_defaultLogger;
    ErrHandler m_errHandler;
    LogLevels m_logLevels;
    std::unique_ptr<details::Formatter> m_formatter;
    Level m_globalLevel = Level::Info;
    Level m_flushLevel = Level::Off;
    bool m_autoRegistration = true;
    size_t m_backtraceN = 0;
};

struct DummyMutex
{
    void lock() const {}
    void unlock() const {}
    bool try_lock() const { return true; }
};

struct SynchronousFactory
{
    template<typename sink_t, typename... Args>
    static std::shared_ptr<Logger> Create(std::string name, Args && ...args)
    {
        auto sink = std::make_shared<sink_t>(std::forward<Args>(args)...);
        auto logger = std::make_shared<Logger>(std::move(name), std::move(sink));
        Registry::Instance().InitializeLogger(logger);
        return logger;
    }

    static std::shared_ptr<Logger> Create(std::string name, SinksInitializerList sinks)
    {
        auto logger = std::make_shared<Logger>(std::move(name), std::move(sinks));
        Registry::Instance().InitializeLogger(logger);
        return logger;
    }
};

//APIs
using FileSink      = sinks::BasicFileSink<DummyMutex>;
using FileSinkMT    = sinks::BasicFileSink<std::mutex>;
using StreamSink    = sinks::BasicOstreamSink<DummyMutex>;
using StreamSinkMT  = sinks::BasicOstreamSink<std::mutex>;

///@brief creates and registers a logger with a templated sink type
template <typename sink_t, typename... Args>
inline std::shared_ptr<Logger> Create(std::string name, Args && ...args)
{
    return SynchronousFactory::Create(std::move(name), std::forward<Args>(args)...);
}

///@brief creates and registers a basic file logger with multi-threading log support
inline std::shared_ptr<Logger> BasicLoggerMT(std::string name, std::string_view filename, bool truncate = false)
{
    return SynchronousFactory::Create<sinks::BasicFileSink<std::mutex> >(std::move(name), filename, truncate);
}

///@brief creates and registers a basic file logger without multi-threading log support
inline std::shared_ptr<Logger> BasicLogger(std::string name, std::string_view filename, bool truncate = false)
{
    return SynchronousFactory::Create<sinks::BasicFileSink<DummyMutex> >(std::move(name), filename, truncate);
}

///@brief creates and registers an out stream logger with multi-threading log support
inline std::shared_ptr<Logger> OstreamLoggerMT(std::string name, std::ostream & os, bool forceFlush = false)
{
    return SynchronousFactory::Create<sinks::BasicOstreamSink<std::mutex> >(std::move(name), os, forceFlush);
}

///@brief creates and registers an out stream logger without multi-threading log support
inline std::shared_ptr<Logger> OstreamLogger(std::string name, std::ostream & os, bool forceFlush = false)
{
    return SynchronousFactory::Create<sinks::BasicOstreamSink<DummyMutex> >(std::move(name), os, forceFlush);
}

///@brief creates and registers a logger with multi sinks
inline std::shared_ptr<Logger> MultiSinksLogger(std::string name, SinksInitializerList sinks)
{
    return SynchronousFactory::Create(std::move(name), std::move(sinks));
}

///@brief initializes a logger
inline void InitializeLogger(std::shared_ptr<Logger> logger)
{
    Registry::Instance().InitializeLogger(std::move(logger));
}

///@brief gets a logger by name
inline std::shared_ptr<Logger> Get(const std::string & name)
{
    return Registry::Instance().Get(name);
}

///@brief gets default logger
inline std::shared_ptr<Logger> DefaultLogger()
{
    return Registry::Instance().DefaultLogger();
}

///@brief gets row pointer of default logger
inline Logger * DefaultLoggerRaw()
{
    return Registry::Instance().GetDefaultRaw();
}

///@brief sets a logger as default logger, note: should not be used concurrently with logging APIs
inline void SetDefaultLogger(std::shared_ptr<Logger> defaultLogger)
{
    Registry::Instance().SetDefaultLogger(std::move(defaultLogger));
}

///@brief dump backtrace messages from default logger
inline void DumpBacktrace()
{
    if (auto log = DefaultLoggerRaw(); log)
        log->DumpBacktrace();
}

///@brief adds trace log message to default log with formatted args
template <typename... Args>
inline void Trace(std::string format, Args && ...args)
{
    if (auto log = DefaultLoggerRaw(); log)
        log->Trace(format, std::forward<Args>(args)...);
}

///@brief adds debug log message to default log with formatted args 
template <typename... Args>
inline void Debug(std::string format, Args && ...args)
{
    if (auto log = DefaultLoggerRaw(); log)
        log->Debug(format, std::forward<Args>(args)...);
}

///@brief adds info log message to default log with formatted args 
template <typename... Args>
inline void Info(std::string format, Args && ...args)
{
    if (auto log = DefaultLoggerRaw(); log)
        log->Info(format, std::forward<Args>(args)...);
}

///@brief adds warn log message to default log with formatted args 
template <typename... Args>
inline void Warn(std::string format, Args && ...args)
{
    if (auto log = DefaultLoggerRaw(); log)
        log->Warn(format, std::forward<Args>(args)...);
}

///@brief adds error log message to default log with formatted args 
template <typename... Args>
inline void Error(std::string format, Args && ...args)
{
    if (auto log = DefaultLoggerRaw(); log)
        log->Error(format, std::forward<Args>(args)...);
}

///@brief adds fator log message to default log with formatted args 
template <typename... Args>
inline void Fatal(std::string format, Args && ...args)
{
    if (auto log = DefaultLoggerRaw(); log)
        log->Fatal(format, std::forward<Args>(args)...);
}

///@brief logs message with specified source and level
inline void Log(details::SourceLoc source, Level level, const std::string & msg)
{
    if (auto log = DefaultLoggerRaw(); log)
        log->Log(source, level, msg);
}

///@brief clears and shuts down all the loggers
inline void ShutDown()
{
    Registry::Instance().ShutDown();
}

}//namespace log
}//namespace generic
