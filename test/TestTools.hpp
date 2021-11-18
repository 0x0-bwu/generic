#ifndef TEST_TESTTOOLSHPP
#define TEST_TESTTOOLSHPP
#define BOOST_TEST_INCLUDED
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include "generic/tools/ProgramOptions.hpp"
#include "generic/tools/Parser.hpp"
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

test_suite * create_tools_test_suite()
{
    test_suite * tools_suite = BOOST_TEST_SUITE("s_tools");
    //
    tools_suite->add(BOOST_TEST_CASE(&t_program_options));
	tools_suite->add(BOOST_TEST_CASE(&t_parser));
    //
    return tools_suite;
}
#endif//TEST_TESTTOOLSHPP