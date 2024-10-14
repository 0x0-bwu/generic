/**
 * @file TestTree.hpp
 * @author bwu
 * @brief Unit test of classes and functions in namespace tree
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/test/TestCommon.hpp"
#include "generic/tree/QuadTreeUtilityMT.hpp"
#include "generic/tree/KdTreeUtilityMT.hpp"
#include "generic/tree/Varification.hpp"
#include "generic/tree/BVHUtilityMT.hpp"
#include "generic/common/Traits.hpp"
#include "generic/tree/IO.hpp"
#include <chrono>
#include <array>
#include <map>
using namespace boost::unit_test;
using namespace generic;
using namespace generic::tree;
using namespace generic::math;
using namespace generic::geometry;
using generic::common::float_type;
using t_tree_num_types = boost::mpl::list<int32_t, int64_t, float, double>;
BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_bvh_collision_t, num_type)
{
    using P = Point3D<float_type<num_type> >;
    using B = Box3D<num_type>;
    using Node = typename bvh::BVH<num_type>::BVHNode;

    struct BoxExtent { B operator() (const B & b) const { return b; } };
    struct BoxCentroid { P operator() (const B & b) const { return b.Center(); } };
    std::vector<B > Objs { {{0, 0, 0}, {2, 2, 2}}, {{1, 1, 1}, {3, 3, 3}}, {{3, 3, 3}, {5, 5, 5}} };

    std::vector<B * > primitives;
    for(size_t i = 0; i < Objs.size(); ++i)
        primitives.push_back(&Objs[i]);
    
    std::vector<B > boxes;
    std::vector<P > centers;
    B boundary = 
    bvh::BVHUtility::CalculateBBoxAndCenter<B, num_type, BoxExtent, BoxCentroid>(primitives, boxes, centers);

    bvh::BVH<num_type> bvh;
    bvh::BinnedSahBuilder<num_type, 16, 1, TopDownTaskSpawner> builder(bvh);
    builder.Build(boundary, boxes, centers);

    std::vector<num_type> coors, t_coors {0,0,0,5,5,5,0,0,0,3,3,3,3,3,3,5,5,5,0,0,0,2,2,2,1,1,1,3,3,3};
    std::vector<size_t> primCount, t_primCount {0, 0, 1, 1, 1};
    std::vector<size_t> firstChilorPrim, t_firstChilorPrim {1, 3, 2, 0, 1};

    for (size_t i = 0; i < bvh.nodeCount; ++i) {
        Node & node = bvh.nodes[i];
        for (size_t j = 0; j < 2; ++j)
            for (size_t k = 0; k < 3; ++k)
                coors.push_back(node.boundary[j][k]);
        
        primCount.push_back(node.primitiveCount);
        firstChilorPrim.push_back(node.firstChildOrPrim);
    }
    BOOST_TEST(t_coors == coors, boost::test_tools::per_element());
    BOOST_TEST(t_primCount == primCount, boost::test_tools::per_element());
    BOOST_TEST(t_firstChilorPrim == firstChilorPrim, boost::test_tools::per_element());

    std::vector<size_t> t_primIndices {0 ,1, 2}; 
    BOOST_TEST(t_primIndices == bvh.primIndices, boost::test_tools::per_element());

    bvh::BVHCollisionDetector<num_type> cd(bvh);
    std::list<std::pair<size_t, size_t> > collisions;
    cd.CollisionDectect(collisions);
    
    std::pair<size_t, size_t> t_collision {0, 1};
    BOOST_TEST(collisions.size() == size_t(1));
    BOOST_CHECK(t_collision == collisions.front());
}

void t_bvh_intersection()
{
    using T = Triangle3D<double>;
    using V = Vector3D<double>;
    using P = Point3D<double>;
    using B = Box3D<double>;

    struct TriExtent { B operator() (const T & t) const { return t.BoundingBox(); } };
    struct TriCentroid { P operator() (const T & t) const { return t.Center(); } };

    std::vector<T> Objs { T(P(0,0,0), P(1,0,1), P(0,1,1)), T(P(0,0,0.5), P(0.5,0,1), P(0, 0.5, 1)) };
    std::vector<T * > primitives;
    for(size_t i = 0; i < Objs.size(); ++i)
        primitives.push_back(&Objs[i]);

    std::vector<B > boxes;
    std::vector<P > centers;
    B boundary = 
    bvh::BVHUtility::CalculateBBoxAndCenter<T, double, TriExtent, TriCentroid>(primitives, boxes, centers);

    bvh::BVH<double> bvh;
    bvh::BinnedSahBuilder<double, 16, 2> builder(bvh);
    builder.Build(boundary, boxes, centers);

    bvh::Ray<double> ray(P(0, 0, 1), V(1, 1, -1));
    bvh::TriangleIntersector<double> triangleIntersector(bvh, primitives);
    bvh::SingleRayTraverser<double> traverser(bvh);
    bvh::SingleRayTraverser<double>::Statistics stat;

    auto hit = traverser.Traverse(ray, triangleIntersector, &stat);
    BOOST_CHECK(hit.Hit(ray));
    BOOST_CHECK(hit.primIndex == 1);
    BOOST_CHECK(hit.intersection.t == 1.0 / 6);
    BOOST_CHECK(hit.intersection.u == 1.0 / 3);
    BOOST_CHECK(hit.intersection.v == 1.0 / 3);
    BOOST_CHECK(stat.intersections == 2);
    BOOST_CHECK(stat.traversalSteps == 0);
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_bvh_t, num_type)
{
    BOOST_TEST_MESSAGE("start bvh test");

    const size_t test_size = 10000;
    using P = Point3D<num_type>;
    using FP = Point3D<float_type<num_type> >;
    using B = Box3D<num_type>;

    struct Extent { B operator() (const B & b) const { return b; } };
    struct Centroid { FP operator() (const B & b) const { return b.Center(); } };
        
    std::vector<B > Objs(test_size);
    auto randP = []{ return P(Random<num_type>(-10000, 10000), Random<num_type>(-10000, 10000), Random<num_type>(-10000, 10000)) ; };
    auto randB = [randP]{ B b(randP(), randP()); b.Normalize(); return b; };
    for (size_t i = 0; i < test_size; ++i) Objs[i] = randB();

    std::vector<B * > primitives(test_size);
    for (size_t i = 0; i < test_size; ++i) {
        primitives[i] = &(Objs[i]);
    }

    BOOST_TEST_MESSAGE("single thread test:");
    {
        auto start = std::chrono::steady_clock::now();
        std::vector<B > boxes;
        std::vector<FP > centers;
        B boundary = 
        bvh::BVHUtility::CalculateBBoxAndCenter<B, num_type, Extent, Centroid>(primitives, boxes, centers);

        bvh::BVH<num_type> bvh;
        bvh::BinnedSahBuilder<num_type, 16, 1, TopDownTaskSpawner> builder(bvh);
        builder.Build(boundary, boxes, centers);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapse = end - start;
        BOOST_TEST_MESSAGE("build time: " + std::to_string(elapse.count()) + "s");

        size_t depth = 0;
        double boxTolerance = 0.0;
        BVHVarification<num_type, B, Extent> varificator(bvh, primitives);
        BOOST_CHECK(varificator.VarifyTreeStructure(boxTolerance, depth) == true);
        BOOST_TEST_MESSAGE("box tolerance: " + std::to_string(boxTolerance / 100) + "%"); 
        BOOST_TEST_MESSAGE("max depth: " + std::to_string(depth)); 
        BOOST_TEST(boxTolerance < 0.011);
    }

    BOOST_TEST_MESSAGE("multi threads test:");
    {
        auto start = std::chrono::steady_clock::now();
        std::vector<B > boxes;
        std::vector<FP > centers;
        B boundary = 
        bvh::BVHUtilityMT::CalculateBBoxAndCenter<B, num_type, Extent, Centroid>(primitives, boxes, centers);

        bvh::BVH<num_type> bvh;
        bvh::BinnedSahBuilder<num_type, 16, 1, TopDownTaskSpawnerMT> builder(bvh);
        builder.Build(boundary, boxes, centers);
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapse = end - start;
        BOOST_TEST_MESSAGE("build time: " + std::to_string(elapse.count()) + "s");

        size_t depth = 0;
        double boxTolerance = 0.0;
        BVHVarification<num_type, B, Extent> varificator(bvh, primitives);
        BOOST_CHECK(varificator.VarifyTreeStructure(boxTolerance, depth) == true);
        BOOST_TEST_MESSAGE("box tolerance: " + std::to_string(boxTolerance / 100) + "%"); 
        BOOST_TEST_MESSAGE("max depth: " + std::to_string(depth));
        BOOST_TEST(boxTolerance < 0.011); 
    }

    BOOST_TEST_MESSAGE("end bvh test");
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_quadtree_t, num_type)
{
    size_t test_sum = 10000;
    struct Extent
    {
        Box2D<num_type> operator()(const Point2D<num_type> & p) const
        {
            return Box2D<num_type>(p, p);
        }
    };

    std::list<Point2D<num_type> > points;
    for (size_t i = 0; i < test_sum; ++i) {
        Point2D<num_type> p(math::Random<num_type>(-10000, 10000), math::Random<num_type>(-10000, 10000));
        points.push_back(p);
    }

    std::list<Point2D<num_type> * > pts1, pts2, out1, out2;
    for(auto & p : points){
        pts1.push_back(&p);
        pts2.push_back(&p);
    }

    using Tree = QuadTree<num_type, Point2D<num_type>, Extent>;
    Tree tree1, tree2;
    tree1.Build(pts1, 1);

    QuadTreeBuilderMT<Point2D<num_type>, Tree> builder(tree2);
    builder.Build(pts2, 1);

    tree1.GetAllObjects(out1);
    tree2.GetAllObjects(out2);

    std::map<int, std::list<Tree * > > allNodes1, allNodes2;
    Tree::GetAllNodesByLevel(&tree1, allNodes1);
    Tree::GetAllNodesByLevel(&tree2, allNodes2);

    BOOST_CHECK(points.size() == out1.size());
    BOOST_CHECK(out1.size() == out2.size());

    auto iter = allNodes1.begin();
    for(; iter != allNodes1.end(); ++iter){
        int level = iter->first;
        BOOST_CHECK(allNodes2.count(level));
        BOOST_CHECK(allNodes1.at(level).size() == allNodes2.at(level).size());
    }
}

BOOST_TEST_CASE_TEMPLATE_FUNCTION(t_kdtree_t, num_type)
{
    BOOST_TEST_MESSAGE("start kdtree test");
    num_type max = std::numeric_limits<num_type>::max();
    max = std::sqrt(max) / 10;
    auto makeVec7 = [max] {
        VectorN<num_type, 7> vec;
        for(size_t i = 0; i < 7; ++i){
            vec[i] = Random(num_type(-max), num_type(max));
        }
        return vec;
    };

    struct Vecterizer
    {
        VectorN<num_type, 7> operator() (const VectorN<num_type, 7> & v) const { return v; }
    };

    size_t total = 1000;
    std::vector<VectorN<num_type, 7> > vectors(total);
    std::vector<VectorN<num_type, 7> * > primitives(total);
    for (size_t i = 0; i < total; ++i) {
        vectors[i] = makeVec7();
        primitives[i] = &vectors[i];
    }

    using namespace kdtree;
    std::array<PlaneSplitMethod, 3> planeSplitMethod{ PlaneSplitMethod::Sequential, PlaneSplitMethod::MaxRange, PlaneSplitMethod::MaxVariance };
    std::array<ValueSplitMethod, 2> valueSplieMethod{ ValueSplitMethod::Random, ValueSplitMethod::Median };

    BOOST_TEST_MESSAGE("single thread test:");
    for (auto plane : planeSplitMethod) {
        for (auto value : valueSplieMethod) {
            BOOST_TEST_MESSAGE("plane split method: " + toString(plane));
            BOOST_TEST_MESSAGE("value split method: " + toString(value));
            auto start = std::chrono::steady_clock::now();
            KdTree<num_type, 7> tree;
            std::vector<VectorN<num_type, 7> > vecs;
            KdTreeUtility::CalculateVector<VectorN<num_type, 7>, 7, num_type, Vecterizer>(primitives, vecs);
            TreeBuilder<num_type, 7, TopDownTaskSpawner> builder(tree, plane, value);
            builder.Build(vecs);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapse = end-start;
            BOOST_TEST_MESSAGE("build time: " + std::to_string(elapse.count()) + "s");

            KdTreeVarification<num_type, 7, VectorN<num_type, 7>, Vecterizer> varificator(tree, primitives);
            size_t depth = 0;
            BOOST_CHECK(varificator.VarifyTreeStructure(depth) == true);
            BOOST_CHECK(varificator.VarifyKNearest(makeVec7(), 7) == true);
            BOOST_CHECK(varificator.VarifyRNearest(makeVec7(), 7) == true);
            BOOST_TEST_MESSAGE("max depth: " + std::to_string(depth));
        }
    }

    BOOST_TEST_MESSAGE("multi threads test:");
    for (auto plane : planeSplitMethod) {
        for (auto value : valueSplieMethod) {
            BOOST_TEST_MESSAGE("plane split method: " + toString(plane));
            BOOST_TEST_MESSAGE("value split method: " + toString(value));
            auto start = std::chrono::steady_clock::now();
            KdTree<num_type, 7> tree;
            std::vector<VectorN<num_type, 7> > vecs;
            KdTreeUtilityMT::CalculateVector<VectorN<num_type, 7>, 7, num_type, Vecterizer>(primitives, vecs);
            TreeBuilder<num_type, 7, TopDownTaskSpawnerMT> builder(tree, plane, value);
            builder.Build(vecs);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapse = end-start;
            BOOST_TEST_MESSAGE("build time: " + std::to_string(elapse.count()) + "s");

            KdTreeVarification<num_type, 7, VectorN<num_type, 7>, Vecterizer> varificator(tree, primitives);
            size_t depth = 0;
            BOOST_CHECK(varificator.VarifyTreeStructure(depth) == true);
            BOOST_CHECK(varificator.VarifyKNearest(makeVec7(), 7) == true);
            BOOST_CHECK(varificator.VarifyRNearest(makeVec7(), 7) == true);
            BOOST_TEST_MESSAGE("max depth: " + std::to_string(depth));
        }
    }
    BOOST_TEST_MESSAGE("end kdtree test");
}

test_suite * create_tree_test_suite()
{
    test_suite * tree_suite = BOOST_TEST_SUITE("s_tree");
    //
    tree_suite->add(BOOST_TEST_CASE_TEMPLATE(t_bvh_collision_t, t_tree_num_types));
    tree_suite->add(BOOST_TEST_CASE(&t_bvh_intersection));
    tree_suite->add(BOOST_TEST_CASE_TEMPLATE(t_bvh_t, t_tree_num_types));
    tree_suite->add(BOOST_TEST_CASE_TEMPLATE(t_kdtree_t, t_tree_num_types));
    tree_suite->add(BOOST_TEST_CASE_TEMPLATE(t_quadtree_t, t_tree_num_types));
    //
    return tree_suite;
}