/**
 * @file TestTools.hpp
 * @author bwu
 * @brief Unit test of classes and functions of tools
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/thread/ThreadPool.hpp"
#include "generic/tools/ProgramOptions.hpp"
#include "generic/tools/Parser.hpp"
#include "generic/tools/Tools.hpp"
using namespace boost::unit_test;
using namespace generic;
void t_program_options()
{
    float f;
	int m, i;
	bool v;

    using namespace generic::program_options;
	OptionParser op("test options");
	auto helpOption =      op.Add<Switch>("h", "help", "produce help message");
	auto boolOption =      op.Add<Switch, Attribute::Optional>("v", "verbose", "be verbose", &v);
	auto hiddenOption =    op.Add<Switch, Attribute::Hidden>("x", "", "hidden option");
	auto doubleOption =    op.Add<Value<double>>("d", "double", "test for double values", 1.12345678910);
	auto floatOption =     op.Add<Value<float>>("f", "float", "test for float values", 2.12345678910f, &f);
	                       op.Add<Value<int>>("i", "int", "test for int value w/o option", 23, &i);
	auto stringOption =    op.Add<Value<std::string>>("s", "string", "test for string values");
	auto impIntOption =    op.Add<Implicit<int>>("m", "implicit", "implicit test", 42);
	auto advancedOption =  op.Add<Switch, Attribute::Advanced>("", "advanced", "advanced option");
	auto expertOption =    op.Add<Switch, Attribute::Expert>("", "expert", "expert option");
	auto inactiveOption =  op.Add<Switch>("", "inactive", "inactive option");
	inactiveOption->SetAttribute(Attribute::Inactive);
	impIntOption->AssignTo(&m);

    //todo
	BOOST_CHECK(true);
}

void t_parser()
{
	using namespace generic::parser;
    //todo
	BOOST_CHECK(true);
}

void t_timer()
{
	using namespace generic::tools;
    
	auto task = []{
        auto t = tools::AccumulatedTimer::InsertTimer();
        tools::SleepMilliseconds(1);
    };

	size_t total{1000};
	tools::Timer timer;
    thread::ThreadPool pool(4);
    for(size_t i = 0; i < total; ++i) {
        pool.Submit(task);
    }
    pool.Wait();
    tools::AccumulatedTimer::SetUnit(unit::Time::Second);
   	auto accumulated = tools::AccumulatedTimer::Count();
	BOOST_CHECK(accumulated.first > timer.Count());
	BOOST_CHECK(accumulated.second == total);
}

test_suite * create_tools_test_suite()
{
    test_suite * tools_suite = BOOST_TEST_SUITE("s_tools");
    //
    tools_suite->add(BOOST_TEST_CASE(&t_program_options));
	tools_suite->add(BOOST_TEST_CASE(&t_parser));
	tools_suite->add(BOOST_TEST_CASE(&t_timer));
    //
    return tools_suite;
}