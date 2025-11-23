/**
 * @file TestMath.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace image
 * @version 0.1
 * @date 2024-09-12
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/tools/FileSystem.hpp"
#include "generic/image/QRCode.hpp"
using namespace boost::unit_test;
using namespace generic;
using namespace generic::img;

void t_image_qrcode()
{
#ifdef GENERIC_BOOST_GIL_IO_PNG_SUPPORT
    using namespace qr;
    std::string s{"hello world"};
    UnsignedBytes raw(s.begin(), s.end());
    std::string filename = GetTestOutDataPath() + "/qrcode.png";
    BOOST_CHECK(EncodeQR(filename, raw, EncodingMode::BYTE, 0, false));
    BOOST_CHECK(fs::FileExists(filename));
#endif //GENERIC_BOOST_GIL_IO_PNG_SUPPORT
}

test_suite * create_image_test_suite()
{
    test_suite * image_suite = BOOST_TEST_SUITE("s_image");
    //
    image_suite->add(BOOST_TEST_CASE(&t_image_qrcode));
    //
    return image_suite;
}
