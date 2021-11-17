#ifndef TEST_TESTTHREAD_HPP
#define TEST_TESTTHREAD_HPP
#define BOOST_TEST_INCLUDED
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include "thread/TaskFlow.hpp"
#include <chrono>
#include <atomic>
using namespace boost::unit_test;
using namespace generic;
using namespace generic::thread;

void t_taskflow()
{
    TaskFlow taskflow;
    std::atomic<int> num(0);

    auto add2 = [&num]{ num.fetch_add(2);};
    auto sub3 = [&num]{ num.fetch_sub(3);};
    auto mul4 = [&num]{ num.exchange(num.load() * 4);};
    auto div5 = [&num]{ num.exchange(num.load() / 5);};

    auto a = taskflow.Emplace(add2);
    auto b = taskflow.Emplace(sub3);
    auto c = taskflow.Emplace(mul4);
    auto d = taskflow.Emplace(div5);
    auto e = taskflow.Submit([&num]{ return num + 43; });

    d->Precede(a, b);
    a->Success(c);
    b->Success(c);
    e.first->Success(a, b);

    Executor executor(4);
    executor.Run(taskflow);
    BOOST_CHECK(e.second.get() == 42);
}

test_suite * create_thread_test_suite()
{
    test_suite * thread_suite = BOOST_TEST_SUITE("s_thread");
    //
    thread_suite->add(BOOST_TEST_CASE(&t_taskflow));
    //
    return thread_suite;
}
#endif//TEST_TESTTHREAD_HPP
