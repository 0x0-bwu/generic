#include "generic/geometry/Utility.hpp"
#include "generic/geometry/GeometryIO.hpp"
void test()
{
    Polygon2D<double> rect;
    rect << Point2D<double>(0, 0) << Point2D<double>(10, 0) << Point2D<double>(10, 10) << Point2D<double>(0, 10);
    auto res = RoundCorners(rect, 2.0);
    std::cout << res << std::endl;
    std::cout << res.Area() << std::endl;
}

int main()
{
    test();
    return 0;
}