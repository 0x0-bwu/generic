/**
 * @file TestGeometry.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace geometry
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/geometry/BoostGeometryRegister.hpp"
#include "generic/geometry/BoostPolygonRegister.hpp"
#include "generic/geometry/BooleanOperation.hpp"
#include "generic/geometry/Connectivity.hpp"
#include "generic/geometry/PolygonMerge.hpp"
#include "generic/geometry/Triangulator.hpp"
#include "generic/geometry/Transform.hpp"
#include "generic/geometry/Clipper.hpp"
#include "generic/geometry/Utility.hpp"
#include "generic/math/Numbers.hpp"
using namespace boost::unit_test;
using namespace generic;
using namespace generic::geometry;
using generic::common::float_type;
using t_geometry_num_types = boost::mpl::list<float, double, int32_t, int64_t>;

void t_geometry_traits()
{
    bool b1 = std::is_same<typename traits::common_coor_t<int, long long>::type, long long>::value;
    bool b2 = std::is_same<typename traits::common_coor_t<double, int, float>::type, float>::value;
    bool b3 = std::is_same<typename traits::common_coor_t<float, double, double, int, short, long>::type, float>::value;
    bool b4 = std::is_same<typename traits::common_coor_t<long double, long long>::type, long double>::value;
    BOOST_CHECK(b1);
    BOOST_CHECK(b2);
    BOOST_CHECK(b3);
    BOOST_CHECK(b4);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_geometry_point_t, num_type)
{
    auto t = 1e-8;
    //point2d
    Point2D<num_type> pt2d;
    BOOST_CHECK(pt2d[0] == 0 && pt2d[1] == 0);

    //operators
    pt2d = Point2D<num_type>(-5, -5);
    pt2d = -pt2d;
    pt2d += Point2D<num_type>(5, 5);
    pt2d -= Point2D<num_type>(5, 5);
    pt2d = pt2d * 100; pt2d *= 0.01;
    pt2d = pt2d / 0.1; pt2d /= 10;
    BOOST_CHECK_CLOSE((float_type<num_type>)pt2d[0], 5, t);
    BOOST_CHECK_CLOSE((float_type<num_type>)pt2d[1], 5, t); 
    //funcs
    BOOST_CHECK_CLOSE((float_type<num_type>)pt2d.Dot(pt2d), 50, t);
    BOOST_CHECK_CLOSE((float_type<num_type>)pt2d.NormSquare(), 50, t);
    BOOST_CHECK_CLOSE((float_type<num_type>)pt2d.CrossProduct(pt2d), 0, t);
  
    //point3d
    Point3D<num_type> pt3d;
    BOOST_CHECK(pt3d[0] == 0 && pt3d[1] == 0 && pt3d[2] == 0);
    //operators
    pt3d = Point3D<num_type>(-5, -5, -5);
    pt3d = -pt3d;
    pt3d += Point3D<num_type>(5, 5, 5);
    pt3d -= Point3D<num_type>(5, 5, 5);
    pt3d = pt3d * 100; pt3d *= 0.01;
    pt3d = pt3d / 0.1; pt3d /= 10;
    BOOST_CHECK_CLOSE((float_type<num_type>)pt3d[0], 5, t);
    BOOST_CHECK_CLOSE((float_type<num_type>)pt3d[1], 5, t);
    BOOST_CHECK_CLOSE((float_type<num_type>)pt3d[2], 5, t);
    //funcs
    BOOST_CHECK_CLOSE((float_type<num_type>)pt3d.Dot(pt3d), 75, t);
    BOOST_CHECK_CLOSE((float_type<num_type>)pt3d.NormSquare(), 75, t);
    BOOST_CHECK_CLOSE((float_type<num_type>)pt3d.CrossProduct(pt3d).NormSquare(), 0, t);
}

void t_geometry_segment()
{
    auto t = 1e-8;
    Point2D<double> p0(-1.0264718499965966, 9.6163341007195407e-7);
    Point2D<double> p1(0.91950808032415809, -1.0094441102690283e-6);
    Point2D<double> q0(-1.0629447383806110, 9.2709540082141753e-7);
    Point2D<double> q1(1.0811583868227901, -1.0670017169567367e-6);

    Segment2D<double> p(p0, p1), q(q0, q1);
    BOOST_CHECK_CLOSE(Segment2D<double>::Distance(p, q), 0.00055025506003679664, t);

    Point3D<double> u0(0.77998990099877119, 0.61192502360790968, -0.22703111823648214);
    Point3D<double> u1(0.53215344529598951, 0.85724585503339767, -0.10102437809109688);
    Point3D<double> v0(-0.21277333982288837, 0.35091548087075353, -0.49557160679250956);
    Point3D<double> v1(0.11881479667499661, 0.022494725417345762, -0.66426620958372951);

    Segment3D<double> u(u0, u1), v(v0, v1);
    BOOST_CHECK_CLOSE(Segment3D<double>::Distance(u, v), 0.98292397116488739, t);

    Point3D<long int> i0(1e9, 1e9, 0), i1(2e9, 2e9, 0);
    Point3D<long int> j0(1e9, 2e9, 0), j1(2e9, 3e9, 0);
    Segment3D<long int> i(i0, i1), j(j0, j1);
    BOOST_CHECK_CLOSE(Segment3D<long int>::Distance(i, j), 5e8 * std::sqrt(2.0), t);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_geometry_triangle_t, num_type)
{
    auto t = 1e-8;
    Triangle2D<num_type> tri2d(Point2D<num_type>(0, 0), Point2D<num_type>(1, 0), Point2D<num_type>(1, 1));
    auto bbox2d = tri2d.BoundingBox();
    BOOST_CHECK_EQUAL(bbox2d[0][0], 0);
    BOOST_CHECK_EQUAL(bbox2d[0][1], 0);
    BOOST_CHECK_EQUAL(bbox2d[1][0], 1);
    BOOST_CHECK_EQUAL(bbox2d[1][1], 1);

    auto center2d = tri2d.Center();
    BOOST_CHECK_CLOSE(center2d[0], float_type<num_type>(2.0 / 3), t);
    BOOST_CHECK_CLOSE(center2d[1], float_type<num_type>(1.0 / 3), t);
    BOOST_CHECK_CLOSE(tri2d.Area(), float_type<num_type>(0.5), t);

    BOOST_CHECK(tri2d.isCCW());
    tri2d[1] = Point2D<num_type>(0 ,1);
    tri2d.Reverse();
    BOOST_CHECK(tri2d.isCCW());

    Triangle3D<num_type> tri3d(Point3D<num_type>(0, 0, 0), Point3D<num_type>(1, 0, 1), Point3D<num_type>(0, 1, 1));
    auto bbox3d = tri3d.BoundingBox();
    BOOST_CHECK_EQUAL(bbox3d[0][0], 0);
    BOOST_CHECK_EQUAL(bbox3d[0][1], 0);
    BOOST_CHECK_EQUAL(bbox3d[0][2], 0);
    BOOST_CHECK_EQUAL(bbox3d[1][1], 1);
    BOOST_CHECK_EQUAL(bbox3d[1][2], 1);
    BOOST_CHECK_EQUAL(bbox3d[1][2], 1);

    auto center3d = tri3d.Center();
    BOOST_CHECK_CLOSE(center3d[0], float_type<num_type>(1.0 / 3), t);
    BOOST_CHECK_CLOSE(center3d[1], float_type<num_type>(1.0 / 3), t);
    BOOST_CHECK_CLOSE(center3d[2], float_type<num_type>(2.0 / 3), t);
    BOOST_CHECK_CLOSE(tri3d.Area(), float_type<num_type>(0.5 * std::sqrt(3.0)), t);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_geometry_box_t, num_type)
{
    auto max = std::numeric_limits<num_type>::max();
    
    //box2d
    Box2D<num_type> box2d_1;
    BOOST_CHECK(!box2d_1.isValid());
    BOOST_CHECK(box2d_1[0][0] == max && box2d_1[1][0] == -max);
    BOOST_CHECK(box2d_1[0][1] == max && box2d_1[1][1] == -max);

    box2d_1[1] = Point2D<num_type>(2, 2);
    Box2D<num_type> box2d_2(1, 1, 3, 3);
    //opeartors
    box2d_1 += box2d_1;
    box2d_1 |= box2d_2;
    box2d_1 |= Point2D<num_type>(-2, -2);
    BOOST_CHECK(box2d_1[0][0] == -2 && box2d_1[1][0] == 3);
    BOOST_CHECK(box2d_1[0][1] == -2 && box2d_1[1][1] == 3);
    BOOST_CHECK(box2d_2 <= box2d_1);
    BOOST_CHECK(box2d_1 >= box2d_2);
    BOOST_CHECK(!(box2d_2 < box2d_1));
    BOOST_CHECK(!(box2d_1 > box2d_2));
    BOOST_CHECK(box2d_1 > Point2D<num_type>()); 

    //funcs
    BOOST_CHECK(box2d_1.Length() == 5);
    BOOST_CHECK(box2d_1.Width() == 5);
    BOOST_CHECK(box2d_1.Area() == 25);
    auto center2d = box2d_1.Center();
    BOOST_CHECK(center2d[0] == 0.5 && center2d[1] == 0.5);
    Box2D<num_type> box2d_3;
    box2d_3.Normalize();
    BOOST_CHECK(box2d_3.isValid());

    //box3d
    Box3D<num_type> box3d_1;
    BOOST_CHECK(!box3d_1.isValid());
    BOOST_CHECK(box3d_1[0][0] == max && box3d_1[1][0] == -max);
    BOOST_CHECK(box3d_1[0][1] == max && box3d_1[1][1] == -max);
    BOOST_CHECK(box3d_1[0][2] == max && box3d_1[1][2] == -max);

    box3d_1[1] = Point3D<num_type>(2, 2, 2);
    Box3D<num_type> box3d_2(1, 1, 1, 3, 3, 3);
    //opeartors
    box3d_1 += box3d_1;
    box3d_1 |= box3d_2;
    box3d_1 |= Point3D<num_type>(-2, -2, -2);
    BOOST_CHECK(box3d_1[0][0] == -2 && box3d_1[1][0] == 3);
    BOOST_CHECK(box3d_1[0][1] == -2 && box3d_1[1][1] == 3);
    BOOST_CHECK(box3d_1[0][2] == -2 && box3d_1[1][2] == 3);
    BOOST_CHECK(box3d_2 <= box3d_1);
    BOOST_CHECK(box3d_1 >= box3d_2);
    BOOST_CHECK(!(box3d_2 < box3d_1));
    BOOST_CHECK(!(box3d_1 > box3d_2));
    BOOST_CHECK(box3d_1 > Point3D<num_type>()); 

    //funcs
    BOOST_CHECK(box3d_1.Length() == 5);
    BOOST_CHECK(box3d_1.Width() == 5);
    BOOST_CHECK(box3d_1.Width() == 5);
    BOOST_CHECK(box3d_1.HalfArea() == 75);
    BOOST_CHECK(box3d_1.SurfArea() == 150);
    BOOST_CHECK(box3d_1.Volume() == 125);
    BOOST_CHECK(box3d_1.LargestAxis() == 0);
    auto diagonal = box3d_1.Diagonal();
    BOOST_CHECK(diagonal[0] == 5 && diagonal[1] == 5 && diagonal[2] == 5);
    auto center3d = box3d_1.Center();
    BOOST_CHECK(center3d[0] == 0.5 && center3d[1] == 0.5 && center3d[2] == 0.5);
    BOOST_CHECK( Box3D<num_type>::Collision(box3d_1, box3d_2));
    BOOST_CHECK(!Box3D<num_type>::Collision(box3d_1, Box3D<num_type>(6, 6, 6, 7, 7, 7)));
    Box3D<num_type> box3d_3;
    box3d_3.Normalize();
    BOOST_CHECK(box3d_3.isValid());
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_geometry_vector_t, num_type)
{
    VectorN<num_type, 5> v1(1);
    VectorN<num_type, 5> v2(2);
    v2.Set(1, 2, 3, 4, 5);
    BOOST_CHECK(v1 != v2);
    BOOST_CHECK(v2[4] == 5);
    BOOST_CHECK_EQUAL(v1.Dot(v2), 15);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_geometry_utility_t, num_type)
{
    auto t = 1e-5;
    //point-plane
    Point3D<num_type> p1(0, 0, 1);
    Plane<num_type> plane(Point3D<num_type>(0, 0, 0), Point3D<num_type>(1, 0, 1), Point3D<num_type>(0, 1, 1));
    auto cp = ClosestPointOnPlane(p1, plane);
    BOOST_CHECK_CLOSE(cp[0], float_type<num_type>(1.0 / 3), t);
    BOOST_CHECK_CLOSE(cp[1], float_type<num_type>(1.0 / 3), t);
    BOOST_CHECK_CLOSE(cp[2], float_type<num_type>(2.0 / 3), t);
    auto d = PointPlaneDistance(p1, plane);
    BOOST_CHECK_CLOSE(d, float_type<num_type>(std::sqrt(3) / 3), t);

    //point-segment
    auto lp3d = ClosestPointOnSegment(Point3D<num_type>(0, 0, 3), Segment3D<num_type>(Point3D<num_type>(0, 0, 0), Point3D<num_type>(3, 3, 3)));
    BOOST_CHECK_CLOSE(lp3d[0], float_type<num_type>(1), t);
    BOOST_CHECK_CLOSE(lp3d[1], float_type<num_type>(1), t);
    BOOST_CHECK_CLOSE(lp3d[2], float_type<num_type>(1), t);

    auto lp2d = ClosestPointOnSegment(Point2D<num_type>(0, 2), Segment2D<num_type>(Point2D<num_type>(0, 0), Point2D<num_type>(2, 2)));
    BOOST_CHECK_CLOSE(lp2d[0], float_type<num_type>(1), t);
    BOOST_CHECK_CLOSE(lp2d[1], float_type<num_type>(1), t);

    auto lpDist3d = PointSegmentDistanceSq(Point3D<num_type>(0, 0, 1), Segment3D<num_type>(Point3D<num_type>(0, 0, 0), Point3D<num_type>(0, 1, 1)));
    auto lpDist2d = PointSegmentDistanceSq(Point2D<num_type>(0, 1), Segment2D<num_type>(Point2D<num_type>(0, 0), Point2D<num_type>(1, 1)));
    BOOST_CHECK_CLOSE(lpDist3d, 0.5, t);
    BOOST_CHECK_CLOSE(lpDist2d, 0.5, t);

    //point-box
    auto bp3d = ClosestPointInBox(Point3D<num_type>(-1, -1, -1), Box3D<num_type>(0, 0, 0, 2, 2, 2));
    BOOST_CHECK_CLOSE(bp3d[0], float_type<num_type>(0), t);
    BOOST_CHECK_CLOSE(bp3d[1], float_type<num_type>(0), t);
    BOOST_CHECK_CLOSE(bp3d[2], float_type<num_type>(0), t);

    auto bp2d = ClosestPointInBox(Point2D<num_type>(-1, -1), Box2D<num_type>(0, 0, 2, 2));
    BOOST_CHECK_CLOSE(bp2d[0], float_type<num_type>(0), t);
    BOOST_CHECK_CLOSE(bp2d[1], float_type<num_type>(0), t);

    auto bpDist3d = PointBoxDistanceSq(Point3D<num_type>(-1, -1, -1), Box3D<num_type>(0, 0, 0, 2, 2, 2));
    auto bpDist2d = PointBoxDistanceSq(Point2D<num_type>(-1, -1), Box2D<num_type>(0, 0, 2, 2));
    BOOST_CHECK_CLOSE(float_type<num_type>(bpDist3d), float_type<num_type>(3), t);
    BOOST_CHECK_CLOSE(float_type<num_type>(bpDist2d), float_type<num_type>(2), t);

    //angle
    Triangle2D<num_type> tri2d(Point2D<num_type>(0, 0), Point2D<num_type>(1, 0), Point2D<num_type>(0, 1));
    BOOST_CHECK_CLOSE(InteriorAngle<0>(tri2d), InteriorAngle(tri2d, 3), t);
    BOOST_CHECK_CLOSE(InteriorAngle<1>(tri2d), InteriorAngle(tri2d, 4), t);
    BOOST_CHECK_CLOSE(InteriorAngle<2>(tri2d), InteriorAngle(tri2d, 5), t);
    BOOST_CHECK_CLOSE(InteriorAngle<6>(tri2d), math::pi_half   , t);
    BOOST_CHECK_CLOSE(InteriorAngle<7>(tri2d), math::pi_quarter, t);
    BOOST_CHECK_CLOSE(InteriorAngle<8>(tri2d), math::pi_quarter, t);

    Triangle3D<num_type> tri3d(Point3D<num_type>(0, 0, 0), Point3D<num_type>(1, 0, 1), Point3D<num_type>(0, 1, 1));
    BOOST_CHECK_CLOSE(InteriorAngle<0>(tri3d), InteriorAngle(tri3d, 3), t);
    BOOST_CHECK_CLOSE(InteriorAngle<1>(tri3d), InteriorAngle(tri3d, 4), t);
    BOOST_CHECK_CLOSE(InteriorAngle<2>(tri3d), InteriorAngle(tri3d, 5), t);
    BOOST_CHECK_CLOSE(InteriorAngle<6>(tri3d), math::pi / 3, t);
    BOOST_CHECK_CLOSE(InteriorAngle<7>(tri3d), math::pi / 3, t);
    BOOST_CHECK_CLOSE(InteriorAngle<8>(tri3d), math::pi / 3, t);

    //intersection
    Segment2D<num_type> s1(Point2D<num_type>(0, 0), Point2D<num_type>(2, 2));
    Line2D<num_type> line1 = makeLineByTwoPoints(Point2D<num_type>(2, 0), Point2D<num_type>(3, -1));
    std::vector<Point2D<num_type> > points;
    Intersection(s1, line1, points);
    BOOST_CHECK(points.size() == 1);
    BOOST_CHECK(points[0] == Point2D<num_type>(1, 1));

    Line2D<num_type> line2 = makeLineByTwoPoints(Point2D<num_type>(-2, -2), Point2D<num_type>(-1, -1));
    Intersection(s1, line2, points);
    BOOST_CHECK(points.size() == 2);
    BOOST_CHECK(points[0] == Point2D<num_type>(0, 0));
    BOOST_CHECK(points[1] == Point2D<num_type>(2, 2));

    {
        Point3D<float_type<num_type> > p;
        Segment3D<num_type> s{{0, 0, 1}, {1, 1, 0}};
        Plane<num_type> plane{{0, 0, 0}, {1, 0, 1}, {0, 1, 1}};
        BOOST_CHECK(Intersection(s, plane, p));
        BOOST_CHECK_CLOSE(p[0], float_type<num_type>(1.0 / 3), t);
        BOOST_CHECK_CLOSE(p[1], float_type<num_type>(1.0 / 3), t);
        BOOST_CHECK_CLOSE(p[2], float_type<num_type>(2.0 / 3), t);
    }

    //predicates
    BOOST_CHECK(PointLineLocation::Left   == GetPointSegmentLocation(Point2D<num_type>(0, 1), Segment2D<num_type>({0, 0}, {3, 3})));
    BOOST_CHECK(PointLineLocation::OnLine == GetPointSegmentLocation(Point2D<num_type>(1, 1), Segment2D<num_type>({0, 0}, {3, 3})));
    BOOST_CHECK(PointLineLocation::Right  == GetPointSegmentLocation(Point2D<num_type>(1, 0), Segment2D<num_type>({0, 0}, {3, 3})));
    BOOST_CHECK(PointLineLocation::Right  == GetPointSegmentLocation(Point2D<num_type>(0, 1), Segment2D<num_type>({3, 3}, {0, 0})));
    BOOST_CHECK(PointLineLocation::OnLine == GetPointSegmentLocation(Point2D<num_type>(1, 1), Segment2D<num_type>({3, 3}, {0, 0})));
    BOOST_CHECK(PointLineLocation::Left   == GetPointSegmentLocation(Point2D<num_type>(1, 0), Segment2D<num_type>({3, 3}, {0, 0})));

    //simplify
    {
        using P = Point2D<num_type>;
        std::vector<Point2D<num_type> > points{ P(822860000, 900000000), P(655310000, 900000000), P(655310000, 442610000), P(655310000, 355610000), P(568310000, 355610000),
                                                P(568310000, 442610000), P(655310000, 442610000), P(655310000, 900000000), P(170000000, 900000000), P(170000000, 100000000), 
                                                P(163301270, 74999999), P(144999999, 56698729),P(120000000,50000000), P(94999999, 56698729), P(76698729, 75000000),
                                                P(70000000, 100000000), P(76698729, 125000000), P(95000000, 143301270), P(120000000, 150000000), P(145000000, 143301270),
                                                P(163301270, 125000000), P(170000000, 100000000), P(170000000, 900000000), P(0, 900000000), P(0, 0),
                                                P(580000000, 0), P(787046666, 241993), P(787046666, 453186666), P(773000000, 453186666), P(691286666, 534898666),
                                                P(690826666, 534898666), P(690853333, 806753333), P(822860000, 806753333), P(822860000, 900000000) };

        Polygon2D<num_type> complex;
        complex.Set(points);

        std::list<Polygon2D<num_type> > holes;
        Simplify(complex, holes);
        BOOST_CHECK(complex.Size() == 13);
        BOOST_CHECK(holes.size() == 2);
    }
}

void t_geometry_utility()
{
    auto t = 1e-8;
    ///2d
    //Cast
    Point2D<double> p2d_double(-0.999, 1.999);
    auto p2d_int = p2d_double.template Cast<int>();
    BOOST_TEST(p2d_int[0] == 0);
    BOOST_TEST(p2d_int[1] == 1);

    //Distance
    auto dist2d = Distance(Point2D<int>(0, 0), Point2D<int>(1, 1));
    BOOST_CHECK_CLOSE(dist2d, std::sqrt(2.0), t);

    //CircumCircle
    auto circle = CircumCircle(Point2D<int>(), Point2D<int>(1, 0), Point2D<int>(1, 1));
    BOOST_CHECK_CLOSE(circle.o[0], 0.5, t);
    BOOST_CHECK_CLOSE(circle.o[1], 0.5, t);
    BOOST_CHECK_CLOSE(circle.r, std::sqrt(0.5), t);

    circle = CircumCircle(Point2D<int>(-1, 0), Point2D<int>(0, 0), Point2D<int>(1, 0));
    BOOST_CHECK_CLOSE(circle.o[0], 0, t);
    BOOST_TEST(circle.o[1] < std::numeric_limits<double>::infinity());
    BOOST_TEST(circle.r < std::numeric_limits<double>::infinity());

    //Inverse
    auto inv_vec2d = Inverse(Vector2D<int>(2, -2));
    BOOST_CHECK_CLOSE(inv_vec2d[0],  0.5, t);
    BOOST_CHECK_CLOSE(inv_vec2d[1], -0.5, t);

    //SafeInverse
    auto safeinv_vec2d = SafeInverse(Vector2D<int>(0, -0));
    BOOST_TEST(safeinv_vec2d[0] <=  std::numeric_limits<double>::infinity());
    BOOST_TEST(safeinv_vec2d[1] >= -std::numeric_limits<double>::infinity());

    //CrossProduct
    auto cp2d = CrossProduct(Vector2D<int>(1, 0), Vector2D<int>(0, 1));
    BOOST_TEST(cp2d == 1);

    ///3d
    //Cast
    Point3D<double> p3d_double(0.001, 1.999, -0.999);
    auto p3d_int = p3d_double.template Cast<int>();
    BOOST_TEST(p3d_int[0] == 0);
    BOOST_TEST(p3d_int[1] == 1);
    BOOST_TEST(p3d_int[2] == 0);

    //Inverse
    auto inv_vec3d = Inverse(Vector3D<int>(2, -2, 0));
    BOOST_CHECK_CLOSE(inv_vec3d[0],  0.5, t);
    BOOST_CHECK_CLOSE(inv_vec3d[1], -0.5, t);
    BOOST_TEST(inv_vec3d[2] == std::numeric_limits<double>::infinity());

    //SafeInverse
    auto safeinv_vec3d = SafeInverse(Vector3D<int>(0, -0, 2));
    BOOST_TEST(safeinv_vec3d[0] <=  std::numeric_limits<double>::infinity());
    BOOST_TEST(safeinv_vec3d[1] >= -std::numeric_limits<double>::infinity());
    BOOST_CHECK_CLOSE(safeinv_vec3d[2], 0.5, t);

    //CrossProduct, DotProduct, Normalize
    Vector3D<int> v3d1(-5, 30, -97), v3d2(62, -42, 0);
    auto cp3d = CrossProduct(v3d1, v3d2);
    BOOST_TEST(DotProduct(v3d1, cp3d) == 0);
    BOOST_TEST(DotProduct(v3d2, cp3d) == 0);

    auto norm = Normalize(cp3d);
    BOOST_CHECK_CLOSE(norm.NormSquare(), 1.0, t);

    //PointLineLocation
    Segment2D<double> seg(Point2D<double>(0.0, 0.0), Point2D<double>(2.0, 2.0));
    BOOST_CHECK(PointLineLocation::OnLine == GetPointSegmentLocation(Point2D<double>(1.0, 1.0), seg));
    BOOST_CHECK(PointLineLocation::Right  == GetPointSegmentLocation(Point2D<double>(2.0, 1.0), seg));
    BOOST_CHECK(PointLineLocation::Left   == GetPointSegmentLocation(Point2D<double>(1.0, 2.0), seg));

    //PointTriangleLocation
    Point2D<int> p1(-2, 0), p2(2, 0), p3(0, 2);
    BOOST_CHECK(PointTriangleLocation::Outside == GetPointTriangleLocation(Point2D<int>( 0, -1), p1, p2, p3));
    BOOST_CHECK(PointTriangleLocation::OnEdge1 == GetPointTriangleLocation(Point2D<int>( 0,  0), p1, p2, p3));
    BOOST_CHECK(PointTriangleLocation::OnEdge2 == GetPointTriangleLocation(Point2D<int>( 1,  1), p1, p2, p3));
    BOOST_CHECK(PointTriangleLocation::OnEdge3 == GetPointTriangleLocation(Point2D<int>(-1,  1), p1, p2, p3));
    BOOST_CHECK(PointTriangleLocation::Inside  == GetPointTriangleLocation(Point2D<int>( 0,  1), p1, p2, p3));
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_triangulation_t, num_type)
{
    using Edge = std::pair<size_t, size_t >;
    using Point = geometry::Point2D<num_type>;
    using Triangulation = geometry::tri::Triangulation<Point>;
    using Triangulator = geometry::tri::Triangulator2D<num_type>;
    using geometry::tri::VertexVec;
    using geometry::tri::TriangleVec;
    std::vector<Point> points { {-4000, 0}, {-3000, -2000}, { 3000, -2000}, { 4000,     0}, { 3000,  2000}, { -3000,  2000}, 
                                {-3000, 0}, {-2000, -1000}, {-1000,     0}, {-2000,  1000}, 
                                { 1000, 0}, { 2000, -1000}, { 3000,     0}, { 2000,  1000}};

    std::vector<Edge > edges { {0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0},
                               {6, 7}, {7, 8}, {8, 9}, {9, 6},
                               {10, 11}, {11,12}, {12, 13}, {13, 10}};

    Triangulation triangulation, triangulation2;
    Triangulator triangulator(triangulation), triangulator2(triangulation2);
    triangulator.InsertVertices(points.begin(), points.end(), [](const Point & p){ return p[0]; }, [](const Point & p){ return p[1]; });
    triangulator.InsertEdges(edges.begin(), edges.end(), [](const Edge & e){ return e.first; }, [](const Edge & e){ return e.second; });
    triangulation2 = triangulation;
    BOOST_CHECK_NO_THROW(triangulator.EraseSuperTriangle());
    VertexVec vertexVec {
                        {0,  {0 ,4}},
                        {1,  {1, 0, 7}},
                        {2,  {12, 2, 1, 15}},
                        {3,  {2, 17}},
                        {4,  {13, 19, 3, 17}},
                        {5,  {3, 10, 4}},
                        {6,  {4, 0, 7, 6, 10, 8}},
                        {7,  {1, 8, 7, 14, 5, 12}},
                        {8,  {5, 9, 6, 8}},
                        {9,  {10, 3, 11, 6, 19, 9}},
                        {10, {5, 9, 16, 11, 18, 14}},
                        {11, {14, 15, 16, 12}},
                        {12, {13, 17, 16, 2, 18, 15}},
                        {13, {13, 18, 11, 19}}};

    auto max = std::numeric_limits<size_t>::max();
    TriangleVec triVec  {
                        {{ 0,  1,  6}, {max,   7,   4}},
                        {{ 1,  2,  7}, {max,  12,   7}},
                        {{ 3, 12,  2}, { 17,  15, max}},
                        {{ 4,  5,  9}, {max,  10,  19}},
                        {{ 0,  6,  5}, {  0,  10, max}},
                        {{ 8,  7, 10}, {  8,  14,   9}},
                        {{ 6,  8,  9}, {  8,   9,  10}},
                        {{ 6,  1,  7}, {  0,   1,   8}},
                        {{ 7,  8,  6}, {  5,   6,   7}},
                        {{ 8, 10,  9}, {  5,  11,   6}},
                        {{ 6,  9,  5}, {  6,   3,   4}},
                        {{ 9, 10, 13}, {  9,  18,  19}},
                        {{ 7,  2, 11}, {  1,  15,  14}},
                        {{12,  4, 13}, { 17,  19,  18}},
                        {{ 7, 11, 10}, { 12,  16,   5}},
                        {{11,  2, 12}, { 12,   2,  16}},
                        {{11, 12, 10}, { 15,  18,  14}},
                        {{ 4, 12,  3}, { 13,   2, max}},
                        {{12, 13, 10}, { 13,  11,  16}},
                        {{ 9, 13,  4}, { 11,  13,   3}}};
    
    for(size_t i = 0; i < triangulation.vertices.size(); ++i){
        //BOOST_CHECK(triangulation.points[triangulation.vertices[i].index] == points[triangulation.vertices[i].index]);
    }
    
    //BOOST_CHECK_EQUAL_COLLECTIONS(triangulation.vertices.begin(), triangulation.vertices.end(), vertexVec.begin(), vertexVec.end());
    
    //BOOST_CHECK_EQUAL_COLLECTIONS(triangulation.triangles.begin(), triangulation.triangles.end(), triVec.begin(), triVec.end());

    // BOOST_CHECK_NO_THROW(triangulator.EraseOuterTriangles());
    BOOST_CHECK_NO_THROW(triangulator2.EraseOuterTrianglesAndHoles());
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_boost_polygon_t, num_type)
{
    auto t = 1e-5;
    using namespace boost::polygon;
    using namespace boost::polygon::operators;

    BOOST_CHECK_CLOSE(euclidean_distance(Point2D<num_type>(0, 0), Point2D<num_type>(1, 1)), std::sqrt(2.0), t);

    Segment2D<num_type> s1(Point2D<num_type>(0, 0), Point2D<num_type>(1, 1));
    Segment2D<num_type> s2(Point2D<num_type>(1, 0), Point2D<num_type>(0, 1));
    BOOST_CHECK(intersects(s1, s2));

    Triangle2D<num_type> tri(Point2D<num_type>(0, 0), Point2D<num_type>(1, 0), Point2D<num_type>(1, 1));
    BOOST_CHECK_CLOSE(area(tri), 0.5, t);

    Box2D<num_type> box(Point2D<num_type>(0, 0), Point2D<num_type>(1, 1));
    BOOST_CHECK_CLOSE(area(box), 1.0, t);

    Polygon2D<num_type> polygon;
    std::vector<Point2D<num_type> > pts { {0, 0}, {1, 0}, {1, 1}, {0, 1} };
    set_points(polygon, pts.begin(), pts.end());
    BOOST_CHECK_CLOSE(perimeter(polygon), 4.0, t);

    auto outline = polygon;
    scale_up(outline, 3);
    std::list<Polygon2D<num_type> > holes {polygon};
    PolygonWithHoles2D<num_type> pwh;
    set_points(pwh, outline.ConstBegin(), outline.ConstEnd());
    set_holes(pwh, holes.begin(), holes.end());
    BOOST_TEST(std::distance(begin_points(pwh), end_points(pwh)) == 4);
    BOOST_TEST(std::distance(begin_holes(pwh), end_holes(pwh)) == 1);
    BOOST_CHECK_CLOSE(area(pwh), 8.0, t);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_boost_geometry_t, num_type)
{
    auto t = 1e-5;
    using namespace boost::geometry;

    BOOST_CHECK_CLOSE(distance(Point2D<num_type>(0, 0), Point2D<num_type>(1, 1)), std::sqrt(2.0), t);

    Segment2D<num_type> s1(Point2D<num_type>(0, 0), Point2D<num_type>(1, 1));
    Segment2D<num_type> s2(Point2D<num_type>(1, 0), Point2D<num_type>(0, 1));
    BOOST_CHECK(intersects(s1, s2));

    Triangle2D<num_type> tri(Point2D<num_type>(0, 0), Point2D<num_type>(1, 0), Point2D<num_type>(1, 1));
    BOOST_CHECK_CLOSE(area(tri), -0.5, t);//pay attention to the direction!

    Box2D<num_type> box(Point2D<num_type>(0, 0), Point2D<num_type>(1, 1));
    BOOST_CHECK_CLOSE(area(box), 1.0, t);

    Polygon2D<num_type> polygon;
    std::vector<Point2D<num_type> > pts { {0, 0}, {0, 3}, {3, 3}, {3, 0} };
    assign_points(polygon, pts);
    BOOST_CHECK_CLOSE(perimeter(polygon), 12.0, t);

    std::vector<Point2D<num_type> > hpts { {1, 1}, {2, 1}, {2, 2}, {1, 2} };
    PolygonWithHoles2D<num_type> pwh;
    assign_points(pwh.outline, pts);
    pwh.holes.resize(1);
    assign_points(pwh.holes.front(), hpts);
    BOOST_CHECK_CLOSE(area(pwh), 8.0, t);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_boolean_operation_t, num_type)
{
    using Point = Point2D<num_type>;
    using Polygon = Polygon2D<num_type>;

    auto t = 1e-5;
    Polygon p1, p2;
    p1 << Point(0,  0) << Point(5, 0) << Point(5, 5) << Point(0, 5);
    p2 << Point(2,  2) << Point(7, 2) << Point(7, 7) << Point(2, 7);

    std::list<Polygon> results;
    boolean::Unite(p1, p2, results);
    BOOST_CHECK(results.size() == 1);
    auto s = boost::polygon::area(results.front());
    BOOST_CHECK_CLOSE(s, 41, t);

    boolean::Subtract(p1, p2, results);
    BOOST_CHECK(results.size() == 1);
    s = boost::polygon::area(results.front());
    BOOST_CHECK_CLOSE(s, 16, t);

    boolean::Intersect(p1, p2, results);
    BOOST_CHECK(results.size() == 1);
    s = boost::polygon::area(results.front());
    BOOST_CHECK_CLOSE(s, 9, t);

    boolean::Xor(p1, p2, results);
    BOOST_CHECK(results.size() == 2);
    auto s1 = boost::polygon::area(results.front());
    auto s2 = boost::polygon::area(results.back());
    BOOST_CHECK_CLOSE(s1 + s2, 32, t);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_connectivity_t, num_type)
{
    using float_t = float_type<num_type>;
    using Point = Point2D<num_type>;
    using PointF = Point2D<float_t>;
    using Triangle = Triangle2D<num_type>;
    using Box = Box2D<num_type>;
    using BoxF = Box2D<float_t>;
    using Segment = Segment2D<num_type>;
    using Polygon = Polygon2D<num_type>;
    using PolygonWithHoles = PolygonWithHoles2D<num_type>;

    using PointSet = std::vector<Point>;
    struct IndexTriangle { std::array<size_t, 3> indices; };
    struct IndexShape { std::vector<size_t> indices; };
    struct IndexInstShape { IndexShape * shape; Point offset; };
    struct Doughnut { Point center; num_type inR, outR; };

    PointSet points {Point(0, 0), Point(30, 0), Point(30, 10), Point(50, 10), Point(50, 20), Point(30, 20)};

    BoxF box{PointF(-2, -2), PointF(2, 2)};
    Doughnut doughnut{Point(0, 0), 20, 30};
    IndexTriangle triangle{std::array<size_t, 3>{0, 1, 2}};
    IndexShape shape{{2, 3, 4, 5}};
    IndexInstShape shape1{&shape, Point(0, 0)}, shape2{&shape, Point(10, 0)}, shape3{&shape, Point(20, 0)};
    std::vector<const IndexInstShape * > shapes{&shape1, &shape2, &shape3};

    auto boxGetter = [](const BoxF * box) { return box->template Cast<num_type>(); };
    auto triGetter = [&points](IndexTriangle & tri)
    {
        return Triangle(points[tri.indices[0]], 
                        points[tri.indices[1]],
                        points[tri.indices[2]]);
    };

    struct DoughnutGeomGetter
    {
        static PolygonWithHoles Get(const Doughnut & doughnut)
        {
            Circle<num_type> in {doughnut.center, doughnut.inR };
            Circle<num_type> out{doughnut.center, doughnut.outR};

            auto iPolygon = InscribedPolygon(in, 12);
            auto oPolygon = CircumscribedPolygon(out, 32);

            PolygonWithHoles pwh;
            pwh.outline = std::move(oPolygon);
            pwh.holes.push_back(std::move(iPolygon));
            return pwh;
        }
    };
    
    std::function<PolygonWithHoles(const Doughnut & doughnut)> doughnutGeomGetter = std::bind(&DoughnutGeomGetter::Get, doughnut);

    struct InstShapeGeomGetter
    {
       explicit InstShapeGeomGetter(PointSet & points) : m_points(points) {}
       Polygon operator() (const IndexInstShape * instShape) const
       {
            Polygon polygon;
            std::vector<Point> ps;
            const auto * shape = instShape->shape;
            for(size_t i : shape->indices) ps.push_back(m_points[i]);
            polygon.Set(std::move(ps));

            auto trans = makeShiftTransform2D<float_t>(instShape->offset);
            return trans * polygon;
       }
    private:
        PointSet & m_points;
    };

    InstShapeGeomGetter instShapeGeomGetter(points);

    ConnectivityExtractor<num_type> extractor;

    //Add Objects
    [[maybe_unused]] auto iBox = extractor.AddObject(0, &box, boxGetter);//0
    [[maybe_unused]] auto iNut = extractor.AddObject(1, doughnut, doughnutGeomGetter);//1
    [[maybe_unused]] auto iTri = extractor.AddObject(2, std::ref(triangle), triGetter);//2
    auto iShapes = extractor.AddObjects(2, shapes.begin(), shapes.end(), instShapeGeomGetter);//3-5
    BOOST_CHECK(iShapes.first  == 3);
    BOOST_CHECK(iShapes.second == 5);

    //Add Layer Connections
    extractor.AddLayerConnection(0, 1);
    extractor.AddLayerConnection(1, 2);
    extractor.AddLayerConnection(2, 3);//non-exist layer will not be extracted.

    size_t threads = 2;
    auto graph = extractor.Extract(threads);
    BOOST_CHECK(graph != nullptr);
    
    std::list<index_t> c1, c2;
    std::vector<std::list<index_t> > cc;
    topology::ConnectedComponent(*graph, iBox, c1);
    topology::ConnectedComponent(*graph, iNut, c2);
    topology::ConnectedComponents(*graph, cc);

    c1.sort(std::less<index_t>());
    c2.sort(std::less<index_t>());
    BOOST_CHECK(!(c1 != std::list<index_t>{0}));
    BOOST_CHECK(!(c2 != std::list<index_t>{1, 2, 3, 4, 5}));

    std::sort(cc.begin(), cc.end(), [](const std::list<index_t> & l1, const std::list<index_t> & l2){ return l1.size() < l2.size(); });
    BOOST_CHECK(!(cc != std::vector<std::list<index_t> >{c1, c2}));

    //jumpwire && null layer test

    extractor.Clear();
    //Add Objects
    iBox = extractor.AddObject(0, &box, boxGetter);//0
    iNut = extractor.AddObject(1, doughnut, doughnutGeomGetter);//1
    auto iJw1 = extractor.AddJumpwire(1, 2, Segment(Point(0, -20), Point(30, 10)));//2
    iShapes = extractor.AddObjects(2, shapes.begin(), shapes.end(), instShapeGeomGetter);//3-5
    BOOST_CHECK(iShapes.first  == 3);
    BOOST_CHECK(iShapes.second == 5);
    auto iJw2 = extractor.AddJumpwire(0, 1, Segment(Point(0, 0),  Point(0, 0)));//6

    auto nullLayer = ConnectivityExtractor<num_type>::nullLayer;
    auto iJw3 = extractor.AddJumpwire(0, nullLayer, Segment(Point(0, 0),  Point(0, 0)));//7
    auto iNull1 = extractor.AddObject(nullLayer, &box, boxGetter);//8
    auto iNull2 = extractor.AddObject(nullLayer, &box, boxGetter);//9

    //Add Layer Connections
    extractor.AddLayerConnection(0, 1);
    extractor.AddLayerConnection(1, 2);
    extractor.AddLayerConnection(2, 3);//non-exist layer will not be extracted.
    extractor.AddLayerConnection(2, nullLayer);//null layer will not be extracted.

    threads = 2;
    graph = extractor.Extract(threads);
    BOOST_CHECK(graph != nullptr);

    std::list<index_t> c3, c4;
    topology::ConnectedComponent(*graph, iNull1, c1);
    topology::ConnectedComponent(*graph, iNull2, c2);
    topology::ConnectedComponent(*graph, iBox, c3);
    topology::ConnectedComponent(*graph, iNut, c4);
    topology::ConnectedComponents(*graph, cc);

    c3.sort(std::less<index_t>());
    c4.sort(std::less<index_t>());
    BOOST_CHECK(!(c1 != std::list<index_t>{8}));
    BOOST_CHECK(!(c2 != std::list<index_t>{9}));
    BOOST_CHECK(!(c3 != std::list<index_t>{0, 6, 7}));
    BOOST_CHECK(!(c4 != std::list<index_t>{1, 2, 3, 4, 5}));
    std::sort(cc.begin(), cc.end(),
            [](const std::list<index_t> & l1, const std::list<index_t> & l2){
                bool res = l1.size() < l2.size();
                if(l1.size() == l2.size() && !l1.empty())
                    res = l1.front() < l2.front();
                return res;
            });
    
    BOOST_CHECK(!(cc != std::vector<std::list<index_t> >{c1, c2, c3, c4}));
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_polygon_merge_t, num_type)
{
    using PorpType = size_t;
    using Settings = typename PolygonMerger<PorpType, num_type>::MergeSettings;
    using PolygonData = typename PolygonMerger<PorpType, num_type>::PolygonData;

    Settings settings;
    settings.kernal = Settings::Kernal::Boost;
    PolygonMerger<PorpType, num_type> merger;
    merger.SetMergeSettings(settings);

    
    Box2D<num_type> box1(Point2D<num_type>(0, 0), Point2D<num_type>(100, 100));
    Box2D<num_type> box2(Point2D<num_type>(20, 20), Point2D<num_type>(80, 80));
    PolygonWithHoles2D<num_type> pwh;
    pwh.outline = toPolygon(box1);
    pwh.holes.emplace_back(toPolygon(box2));

    Box2D<num_type> box3(Point2D<num_type>(40, 0), Point2D<num_type>(60, 100));

    merger.AddObject(0, pwh);
    merger.AddObject(0, box3);

    merger.Merge();

    std::list<PolygonData * > polygons;
    merger.GetAllPolygons(polygons);
    BOOST_CHECK(polygons.size() == 1);
    BOOST_CHECK(polygons.front()->holes.size() == 2);

    //parallel merge
    merger.Clear();
    merger.AddObject(0, pwh);
    merger.AddObject(0, box3);

    PolygonMergeRunner runner(merger, 4);
    runner.Run();
    
    polygons.clear();
    merger.GetAllPolygons(polygons);
    BOOST_CHECK(polygons.size() == 1);
    BOOST_CHECK(polygons.front()->holes.size() == 2);

    merger.Clear();
    merger.AddObject(0, pwh);
    merger.AddObject(0, box3);

    //use clipper kernal
    settings.kernal = Settings::Kernal::Clipper;
    merger.SetMergeSettings(settings);
    merger.Merge();

    polygons.clear();
    merger.GetAllPolygons(polygons);
    BOOST_CHECK(polygons.size() == 1);
    BOOST_CHECK(polygons.front()->holes.size() == 2);
}

void t_geometry_additional()
{
    Point2D<double> p1(0, 0), p2(3e-16, -3e-16), p3(-1e-16, 1e-16);
    BOOST_CHECK(p1 != p2);
    BOOST_CHECK(p1 == p3);

    Point3D<double> p4(0, 0, 0), p5(0, 0, -3e-16), p6(-1e-16, 1e-16, 0);
    BOOST_CHECK(p4 != p5);
    BOOST_CHECK(p4 == p6);
}

test_suite * create_geometry_test_suite()
{
    test_suite * geometry_suite = BOOST_TEST_SUITE("s_geometry");
    //
    geometry_suite->add(BOOST_TEST_CASE(&t_geometry_traits));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_geometry_point_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE(&t_geometry_segment));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_geometry_triangle_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_geometry_box_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_geometry_vector_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_geometry_utility_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE(&t_geometry_utility));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_triangulation_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_boost_polygon_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_boost_geometry_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_boolean_operation_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_connectivity_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE_TEMPLATE(t_polygon_merge_t, t_geometry_num_types));
    geometry_suite->add(BOOST_TEST_CASE(&t_geometry_additional));
    //
    return geometry_suite;
}