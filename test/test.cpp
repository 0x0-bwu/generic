#include "generic/image/QRCode.hpp"
#include <iostream>
using namespace generic::img;

int main()
{
    using namespace qr;
    std::string s{"hello world"};
    UnsignedBytes raw(s.begin(), s.end());
    EncodeQR("./test.png", raw, EncodingMode::BYTE, 0, false);
}
