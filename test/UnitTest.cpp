/**
 * @file UnitTest.cpp
 * @author bwu
 * @brief Main entrance of unit test
 * @version 0.1
 * @date 2022-02-22
 */
#define GENERIC_UNIT_TEST
#include "TestMath.hpp"
#include "TestThread.hpp"
#include "TestGeometry.hpp"
#include "TestTree.hpp"
#include "TestTools.hpp"
#include "TestCircuit.hpp"
#include "TestImage.hpp"
#include "TestBoolean.hpp"
#include "TestGraph.hpp"
#include "TestUtils.hpp"
using namespace boost::unit_test;

extern test_suite * create_math_test_suite();
extern test_suite * create_geometry_test_suite();
extern test_suite * create_tree_test_suite();
extern test_suite * create_thread_test_suite();
extern test_suite * create_tools_test_suite();
extern test_suite * create_circuit_test_suite();
extern test_suite * create_image_test_suite();
extern test_suite * create_boolean_test_suite();
extern test_suite * create_graph_test_suite();
extern test_suite * create_utils_test_suite();

void t_additional()
{
    //add additional test here
    BOOST_CHECK(true);
}

test_suite *
init_unit_test_suite(int argc, char* argv[])
{
    framework::master_test_suite().add(create_math_test_suite());
    framework::master_test_suite().add(create_geometry_test_suite());
    framework::master_test_suite().add(create_tree_test_suite());
    framework::master_test_suite().add(create_thread_test_suite());
    framework::master_test_suite().add(create_tools_test_suite());
    framework::master_test_suite().add(create_circuit_test_suite());
    framework::master_test_suite().add(create_image_test_suite());
    framework::master_test_suite().add(create_boolean_test_suite());
    framework::master_test_suite().add(create_graph_test_suite());
    framework::master_test_suite().add(create_utils_test_suite());
    framework::master_test_suite().add(BOOST_TEST_CASE(&t_additional));
    return 0;
}