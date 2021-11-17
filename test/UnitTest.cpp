#define GENERIC_UNIT_TEST
#include <boost/test/included/unit_test.hpp>
#include "TestMath.hpp"
#include "TestThread.hpp"
#include "TestGeometry.hpp"
#include "TestTree.hpp"
#include "TestTools.hpp"
using namespace boost::unit_test;

extern test_suite * create_math_test_suite();
extern test_suite * create_geometry_test_suite();
extern test_suite * create_tree_test_suite();
extern test_suite * create_thread_test_suite();
extern test_suite * create_tools_test_suite();

void t_additional()
{
    //add additional test here
    BOOST_CHECK(true);
}

test_suite *
init_unit_test_suite( int argc, char* argv[] )
{
    framework::master_test_suite().add(create_math_test_suite());
    framework::master_test_suite().add(create_geometry_test_suite());
    framework::master_test_suite().add(create_tree_test_suite());
    framework::master_test_suite().add(create_thread_test_suite());
    framework::master_test_suite().add(BOOST_TEST_CASE(&t_additional));
    return 0;
}