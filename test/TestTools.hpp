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
#include "generic/tools/StringHelper.hpp"
#include "generic/tools/Parser.hpp"
#include "generic/tools/Tools.hpp"
#include "generic/tools/Hash.hpp"
#include "generic/tools/Color.hpp"
#include "generic/tools/Units.hpp"
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

void t_string_helper()
{
	// unit test of WildcardMatch
	BOOST_CHECK(generic::str::WildcardMatch("hello world", "hello*"));
	BOOST_CHECK(generic::str::WildcardMatch("hello world", "*world*"));
	// unit tet of WildcardMatch with ? and *
	BOOST_CHECK(generic::str::WildcardMatch("hello world", "hello?world"));
	BOOST_CHECK(generic::str::WildcardMatch("hello world", "hello?world*"));
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
    tools::ThreadTimer::SetUnit(unit::Time::SECOND);
    tools::AccumulatedTimer<Tag>::SetUnit(unit::Time::SECOND);
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
	BOOST_CHECK(HashCombine(0, "hello world", 42, math::pi) == HashCombine(0, "hello world", 42, math::pi));
	BOOST_CHECK(HashCombine(0, "hello world", 42, math::pi) !=
				HashCombine(0, "hello world", 42, math::pi + std::numeric_limits<float>::epsilon()));	 

	std::vector<double> v1(100);
	std::for_each(v1.begin(), v1.end(), [](auto & d) { d = math::Random<double>(-1e3, 1e3); });
	std::vector<double> v2 = v1;
	std::shuffle(v2.begin(), v2.end(), std::mt19937{std::random_device{}()});
	std::vector<double> v3 = v1; v3[1] += 1e-3;
    BOOST_CHECK(OrderedHash(v1) != OrderedHash(v2));
	BOOST_CHECK(OrderedHash(v1) != OrderedHash(v3));
	BOOST_CHECK(OrderedHash(v2) != OrderedHash(v3));
	BOOST_CHECK(UnorderedHash(v1) == UnorderedHash(v2));
	BOOST_CHECK(UnorderedHash(v1) != UnorderedHash(v3));
}

void t_color()
{
	using namespace generic::color;
	
	// Test predefined colors
	BOOST_CHECK_EQUAL(black, 0x00000000);
	BOOST_CHECK_EQUAL(white, 0xFFFFFFFF);
	BOOST_CHECK_EQUAL(red, 0xFFFF0000);
	BOOST_CHECK_EQUAL(green, 0xFF00FF00);
	BOOST_CHECK_EQUAL(blue, 0xFF0000FF);
	
	// Test RGBToInt
	int color1 = RGBToInt(255, 128, 64);
	BOOST_CHECK((color1 & 0x00FF0000) >> 16 == 255);
	BOOST_CHECK((color1 & 0x0000FF00) >> 8 == 128);
	BOOST_CHECK((color1 & 0x000000FF) == 64);
	
	// Test RGBaToInt
	int color2 = RGBaToInt(255, 128, 64, 200);
	BOOST_CHECK(((color2 >> 24) & 0xFF) == 200);
	
	// Test RGBFromInt
	int r, g, b;
	RGBFromInt(color1, r, g, b);
	BOOST_CHECK_EQUAL(r, 255);
	BOOST_CHECK_EQUAL(g, 128);
	BOOST_CHECK_EQUAL(b, 64);
	
	// Test RGBaFromInt
	int a;
	RGBaFromInt(color2, r, g, b, a);
	BOOST_CHECK_EQUAL(r, 255);
	BOOST_CHECK_EQUAL(g, 128);
	BOOST_CHECK_EQUAL(b, 64);
	BOOST_CHECK_EQUAL(a, 200);
	
	// Test RandomRGB (just ensure it doesn't crash and produces valid values)
	RandomRGB(r, g, b);
	BOOST_CHECK(r >= 0 && r <= 255);
	BOOST_CHECK(g >= 0 && g <= 255);
	BOOST_CHECK(b >= 0 && b <= 255);
	
	// Test toJet
	toJet(0.0, r, g, b);
	BOOST_CHECK_EQUAL(r, 0);
	BOOST_CHECK_EQUAL(g, 0);
	BOOST_CHECK(b >= 127 && b <= 128);
	
	toJet(0.5, r, g, b);
	BOOST_CHECK(r >= 0 && r <= 255);
	BOOST_CHECK(g >= 0 && g <= 255);
	BOOST_CHECK(b >= 0 && b <= 255);
	
	toJet(1.0, r, g, b);
	BOOST_CHECK(r >= 0 && r <= 255);
	BOOST_CHECK(g >= 0 && g <= 255);
	BOOST_CHECK(b >= 0 && b <= 255);
}

void t_units()
{
	using namespace generic::unit;
	
	// Test Time scale conversions
	BOOST_CHECK_CLOSE(Scale2Second(Time::SECOND), 1.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Second(Time::MILLISECOND), 0.001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Second(Time::MICROSECOND), 0.000001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Second(Time::NANOSECOND), 0.000000001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Second(Time::MINUTE), 60.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Second(Time::HOUR), 3600.0f, 0.0001f);
	
	BOOST_CHECK_CLOSE(Scale2Millisecond(Time::SECOND), 1000.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Millisecond(Time::MILLISECOND), 1.0f, 0.0001f);
	
	// Test Length scale conversions
	BOOST_CHECK_CLOSE(Scale2Meter(Length::METER), 1.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Meter(Length::MILLIMETER), 0.001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Meter(Length::MICROMETER), 0.000001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Meter(Length::NANOMETER), 0.000000001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2Meter(Length::INCH), 0.0254f, 0.0001f);
	
	// Test Temperature conversions
	BOOST_CHECK_CLOSE(Celsius2Kelvins(0.0f), 273.15f, 0.0001f);
	BOOST_CHECK_CLOSE(Celsius2Kelvins(100.0f), 373.15f, 0.0001f);
	BOOST_CHECK_CLOSE(Kelvins2Celsius(273.15f), 0.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Kelvins2Celsius(373.15f), 100.0f, 0.0001f);
	
	// Test Scale2SI for various units
	BOOST_CHECK_CLOSE(Scale2SI(Capacitance::F), 1.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Capacitance::PF), 1e-12f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Resistance::OHM), 1.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Resistance::KOHM), 1000.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Current::A), 1.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Current::MA), 0.001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Voltage::V), 1.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Voltage::MV), 0.001f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Power::W), 1.0f, 0.0001f);
	BOOST_CHECK_CLOSE(Scale2SI(Power::MW), 0.001f, 0.0001f);
	
	// Test toString for Time
	BOOST_CHECK_EQUAL(generic::toString(Time::PICOSECOND), "ps");
	BOOST_CHECK_EQUAL(generic::toString(Time::NANOSECOND), "ns");
	BOOST_CHECK_EQUAL(generic::toString(Time::MICROSECOND), "us");
	BOOST_CHECK_EQUAL(generic::toString(Time::MILLISECOND), "ms");
	BOOST_CHECK_EQUAL(generic::toString(Time::SECOND), "sec");
	BOOST_CHECK_EQUAL(generic::toString(Time::MINUTE), "min");
	BOOST_CHECK_EQUAL(generic::toString(Time::HOUR), "hour");
}

test_suite * create_tools_test_suite()
{
    test_suite * tools_suite = BOOST_TEST_SUITE("s_tools");
    //
    tools_suite->add(BOOST_TEST_CASE(&t_program_options));
	tools_suite->add(BOOST_TEST_CASE(&t_string_helper));
	tools_suite->add(BOOST_TEST_CASE(&t_parser));
	tools_suite->add(BOOST_TEST_CASE(&t_timer));
	tools_suite->add(BOOST_TEST_CASE(&t_hash));
	tools_suite->add(BOOST_TEST_CASE(&t_color));
	tools_suite->add(BOOST_TEST_CASE(&t_units));
    //
    return tools_suite;
}