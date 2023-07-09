/**
 * @file TestThread.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace thread
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/thread/MapReduce.hpp"
#include "generic/thread/TaskFlow.hpp"
#include <sstream>
using namespace boost::unit_test;
using namespace generic;
using namespace generic::thread;

namespace primcalc {

using namespace mapreduce;

struct PrimCalcMapTask : public MapTask<long, std::pair<long, long> >
{
    static bool isPrim(long number)
    {
        if(number > 2){
            if(number % 2 == 0) return false;
            long n = std::abs(number);
            long sqrtN = static_cast<long>(std::sqrt(static_cast<double>(n)));
            for(long i = 3; i < sqrtN; i += 2){
                if(n % i == 0) return false;
            }
        }
        else if(number == 0 || number == 1) return false;
        return true;
    }

    template <typename Runner>
    void operator() (Runner & runner, const Key &, const Value & value) const
    {
        for(auto i = value.first; i <= value.second; ++i)
            runner.EmitIntermediate(isPrim(i), i);
    }
};

struct PrimCalcReduceTask : public ReduceTask<bool, long>
{
    template <typename Runner, typename Iterator>
    void operator() (Runner & runner, const Key & key, Iterator begin, Iterator end) const
    {
        if(key)
            std::for_each(begin, end, std::bind(&Runner::Emit, &runner, true, std::placeholders::_1));
    }
};

struct NumberSource
{
    long sequence, step;
    long first, last;
    NumberSource(long _first, long _last, long _step)
        : sequence(0), step(_step), first(_first), last(_last)
    {
    }

    bool SetupKey(long & k)
    {
        k = sequence++;
        return k * step <= last;
    }

    bool GetData(const long & key, std::pair<long, long> & value)
    {
        std::pair<long, long> tmp;
        tmp.first = first + (key * step);
        tmp.second = std::min(tmp.first + step - 1, last);
        
        std::swap(tmp, value);
        return true;
    }
};

}//namespace primcalc

void t_mapreduce()
{
    using namespace primcalc;
    using namespace mapreduce;
    using PrimCalcJob = Job<PrimCalcMapTask, PrimCalcReduceTask, NumberSource>;

    Specification spec;
    spec.mapTasks = 5;
    spec.reduceTasks = 5;

    long primLimit = 1e5;
    NumberSource source(0, primLimit, primLimit / spec.reduceTasks);
    PrimCalcJob job(source, spec);
    Results results;

    // job.Run<schedule::Sequential<PrimCalcJob> >(results);
    job.Run<schedule::Parallel<PrimCalcJob> >(results);

    std::stringstream ss;
    ss << "mapreduce test:"  << GENERIC_DEFAULT_EOL;
    ss << "total run time: " << results.jobRuntime.count()        << "s" << GENERIC_DEFAULT_EOL;
    ss << "mapped tasks: "   << results.counters.mapKeysCompleted << "/" << results.counters.mapKeysExecuted << GENERIC_DEFAULT_EOL;
    ss << "map run time: "   << results.mapRuntime.count()        << "s" << GENERIC_DEFAULT_EOL;
    for(size_t i = 0; i < results.mapTimes.size(); ++i){
        ss << "No." << i + 1 << ", time: " << results.mapTimes.at(i).count() << "s" << GENERIC_DEFAULT_EOL;
    }

    ss << "reduced tasks: "   << results.counters.reduceKeysCompleted << "/" << results.counters.reduceKeysExecuted << GENERIC_DEFAULT_EOL;
    ss << "reduce run time: " << results.reduceRuntime.count()        << "s" << GENERIC_DEFAULT_EOL;
    for(size_t i = 0; i < results.reduceTimes.size(); ++i){
        ss << "No." << i + 1 << ", time: " << results.reduceTimes.at(i).count() << "s" << GENERIC_DEFAULT_EOL;
    }
    BOOST_TEST_MESSAGE(ss.str());
    BOOST_CHECK(true);
}

void t_taskflow()
{
    using namespace taskflow;
    TaskFlow taskflow;
    std::atomic<int> num(0);

    auto add2 = [&num]{ num.fetch_add(2);};
    auto sub3 = [&num]{ num.fetch_sub(3);};
    auto mul4 = [&num]{ num.exchange(num.load() * 4);};
    auto div5 = [&num]{ num.exchange(num.load() / 5);};

    auto a = taskflow.Emplace(add2, "add2");
    auto b = taskflow.Emplace(sub3, "sub3");
    auto c = taskflow.Emplace(mul4, "mul4");
    auto d = taskflow.Emplace(div5, "div5");
    auto e = taskflow.Submit([&num]{ return num + 43; }, "+43");

    d->Precede(a, b);
    a->Success(c);
    b->Success(c);
    e.first->Success(a, b);

    Executor executor(4);
    auto res = executor.Run(taskflow);
    BOOST_CHECK(e.second.get() == 42);

    std::stringstream ss;
    taskflow.PrintTaskGraph(ss);
    BOOST_TEST_MESSAGE("task graph:\n" + ss.str());

    BOOST_CHECK(res);
    res = executor.Run(taskflow);
    BOOST_CHECK(res == false);
}

test_suite * create_thread_test_suite()
{
    test_suite * thread_suite = BOOST_TEST_SUITE("s_thread");
    //
    thread_suite->add(BOOST_TEST_CASE(&t_mapreduce));
    thread_suite->add(BOOST_TEST_CASE(&t_taskflow));
    //
    return thread_suite;
}