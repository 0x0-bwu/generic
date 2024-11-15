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
#include "generic/tools/Hash.hpp"
#include "generic/math/Numbers.hpp"
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

struct Tag { inline static constexpr size_t groups = 2; };

void t_timer()
{
	using namespace generic::tools;
    
	auto task = []{
        auto t1 = tools::ThreadTimer::InsertTimer();
        tools::SleepMilliseconds(1);
        auto t2 = tools::AccumulatedTimer<Tag>::InsertTimer(1);
        tools::SleepMilliseconds(1);
    };

	size_t total{1000};
	tools::Timer timer;
    thread::ThreadPool pool(4);
    for(size_t i = 0; i < total; ++i)
        pool.Submit(task);

    pool.Wait();
    tools::ThreadTimer::SetUnit(unit::Time::Second);
    tools::AccumulatedTimer<Tag>::SetUnit(unit::Time::Second);
   	auto accumulated1 = tools::ThreadTimer::Count(0);
	auto accumulated2 = tools::AccumulatedTimer<Tag>::Count(1);
	BOOST_CHECK(accumulated1.first > accumulated2.first);
	BOOST_CHECK(accumulated2.first > timer.Count());
	BOOST_CHECK(accumulated1.second == total);
	BOOST_CHECK(accumulated2.second == total);
}

void t_hash()
{
	using namespace generic;
	using namespace generic::hash;
	BOOST_CHECK(HashCombine("hello world", 42, math::pi) == HashCombine("hello world", 42, math::pi));
	BOOST_CHECK(HashCombine("hello world", 42, math::pi) !=
				HashCombine("hello world", 42, math::pi + std::numeric_limits<float>::epsilon()));	 
}

test_suite * create_tools_test_suite()
{
    test_suite * tools_suite = BOOST_TEST_SUITE("s_tools");
    //
    tools_suite->add(BOOST_TEST_CASE(&t_program_options));
	tools_suite->add(BOOST_TEST_CASE(&t_parser));
	tools_suite->add(BOOST_TEST_CASE(&t_timer));
	tools_suite->add(BOOST_TEST_CASE(&t_hash));
    //
    return tools_suite;
}