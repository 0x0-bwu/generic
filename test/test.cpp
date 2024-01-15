#include "generic/geometry/Utility.hpp"
#include "generic/geometry/GeometryIO.hpp"

using namespace generic::geometry;

void test()
{
    Point2D<double> fp1, fp2;
    Point2D<int64_t> p0(5, -5), p1(0, 0), p2(3, 3), p3(8, 8);
    std::cout << InscribedCircle<int64_t>(p0, p1, p2, 1.0, fp1, fp2) << std::endl;
    std::cout << "fp1: " << fp1 << ", fp2: " << fp2 << std::endl;
    std::cout << InscribedCircle<int64_t>(p0, p1, p3, 1.0, fp1, fp2) << std::endl;
    std::cout << "fp1: " << fp1 << ", fp2: " << fp2 << std::endl;
}

int main()
{
    test();
    return 0;
}