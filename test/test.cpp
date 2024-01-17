#include "generic/geometry/Utility.hpp"
#define BOOST_GIL_IO_PNG_SUPPORT 1
#include "generic/geometry/GeometryIO.hpp"
void test()
{
    Polygon2D<double> rect;
    std::vector<Point2D<double> > points{{-14.2, -10.4}, {15.45, -10.4}, {15.45, -9.6}, {13.95, -9.6},
        {13.35, -10}, {-13.8, -10}, {-13.8, -2.55}, {-11.1, -2.55}, {-11.1, 11.45}, {-11.85, 11.45}, {-11.85, -2}, {-14.2, -2}};
    rect.Set(std::move(points));
    std::vector<Polygon2D<double> > polygons{rect};
    GeometryIO::WritePNG("./origin.png", polygons.begin(), polygons.end());
    std::vector<Polygon2D<double>> polygons2{RoundCorners(rect, 0.25)};
    GeometryIO::WritePNG("./round.png", polygons2.begin(), polygons2.end());
}

int main()
{
    test();
    return 0;
}