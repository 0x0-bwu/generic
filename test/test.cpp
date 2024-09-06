#include "generic/image/QRCode.hpp"
#include <iostream>
using namespace generic::img;

int main()
{
    using namespace qr;
    std::string s{"hello world"};
    const UnsignedByte* raw = (const UnsignedByte*)s.c_str();
    EncodeQR("./test.png", raw, EncodingMode::BYTE, 0, false);
}