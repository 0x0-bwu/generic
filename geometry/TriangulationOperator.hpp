/**
 * @file TriangulationOperator.hpp
 * @author bwu
 * @brief Triangulation operator that manipulating the triangulation data
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include <boost/geometry/index/rtree.hpp>
#include <boost/container/flat_set.hpp>
#include "Triangulation.hpp"
#include <unordered_map>
#include <unordered_set>
namespace generic  {
namespace geometry {
namespace tri {

///@base template of vertex index tree
namespace bgi = boost::geometry::index;
template <typename point_t>
class VertexIndexTree
{
};

///@brief represents a class that manage/search Point2D
template <typename num_type>
class VertexIndexTree<Point2D<num_type> >
{
    using Point = Point2D<num_type>;
    using IndexPoint = std::pair<Point, VerIdx>;
    bgi::rtree<IndexPoint, bgi::linear<16, 4> > m_rtree;
public:
    void Clear() { m_rtree.clear(); }
    void AddOneVertex(const Point & pos, const VerIdx index)
    {
        m_rtree.insert(std::make_pair(pos, index));
    }

    void RemoveVertex(const Point & pos, const VerIdx index)
    {
        m_rtree.remove(std::make_pair(pos, index));
    }

    VerIdx NearestVertex(const Point & pos) const
    {
        std::vector<IndexPoint> query;
        m_rtree.query(bgi::nearest(pos, 1), std::back_inserter(query));
        return query.front().second;
    }
};

///@brief represents a class that operate triangulation data
template <typename point_t>
class TriangulationOperator
{
    using Point = point_t;
    using Triangle = triangle_type<point_t>;
    using VerIdxSet = std::unordered_set<VerIdx>;
    using TriIdxSet = std::unordered_set<TriIdx>;
    using VerIdxMap = std::unordered_map<VerIdx, VerIdx>;
    using TriIdxMap = std::unordered_map<TriIdx, TriIdx>;
    using TriIdxFlatSet = boost::container::flat_set<TriIdx>;
    using EdgeSet = typename Triangulation<Point>::EdgeSet;
    using Utility = TriangulationUtility<Point>;

    VerIdxSet m_removedVertices;
    TriIdxSet m_removedTriangles;
    Triangulation<Point> & m_triangulation;
    VertexIndexTree<Point> m_vertexIndexTree;
public:
    ///@brief constructs a TriangulationOperator that manipulating the triangulation data
    TriangulationOperator(Triangulation<Point> & triangulation)
     : m_triangulation(triangulation)
    {
        RebuildVertexIndexTree();
    }

    bool isVertexRemoved(VerIdx iv) const { return m_removedVertices.count(iv); }
    bool isTriangleRemoved(TriIdx it) const { return m_removedTriangles.count(it); }
    size_t CurrentTriangleSize() const { return m_triangulation.triangles.size() - m_removedTriangles.size(); }
    size_t CurrentVertexSize() const { return m_triangulation.vertices.size() - m_removedVertices.size(); }
    const TriIdxSet & GetRemovedTriangles() const { return m_removedTriangles; }
    const VerIdxSet & GetRemovedVertices() const { return m_removedVertices; }
    VerIdx AddOneVertex(const Point & p, TriIdxSet triangles = TriIdxSet{})
    {
        if(m_removedVertices.size()){
            VerIdx iv = *m_removedVertices.begin();
            m_removedVertices.erase(iv);
            auto & vertex = m_triangulation.vertices[iv];
            vertex.triangles = triangles;
            m_triangulation.points[vertex.index] = p;
            m_vertexIndexTree.AddOneVertex(p, iv);
            return iv;
        }

        IndexVertex vertex;
        vertex.index = m_triangulation.points.size();
        m_triangulation.points.push_back(p);
        vertex.triangles = triangles;
        VerIdx iv = m_triangulation.vertices.size();
        m_triangulation.vertices.push_back(vertex);
        m_vertexIndexTree.AddOneVertex(p, iv);
        return iv;
    }

    void RemoveOneVertex(VerIdx iv)
    {
        m_removedVertices.insert(iv);
        m_triangulation.vertices[iv].Clear();
        m_vertexIndexTree.RemoveVertex(m_triangulation.VerIdxPoint(iv), iv);
    }

    TriIdx AddOneTriangle()
    {
        TriIdx it;
        if(m_removedTriangles.size()){
            it = *m_removedTriangles.begin();
            m_removedTriangles.erase(it);
        }
        else {
            it = m_triangulation.triangles.size();
            m_triangulation.triangles.push_back(IndexTriangle()); 
        }
        return it;
    }

    TriIdx AddOneTriangle(IndexTriangle t)
    {
        TriIdx it = AddOneTriangle();
        m_triangulation.triangles[it] = t;
        return it;
    }

    template <typename std::enable_if<traits::is_2d_point_t<Point>::value, bool>::type = true>
    TriIdx AddOneTriangle(VerIdx iv1, VerIdx iv2, VerIdx iv3)
    {
        if(!Utility::isCCW(m_triangulation, iv1, iv2, iv3))
            std::swap(iv2, iv3);

        TriIdx it = AddOneTriangle();

        m_triangulation.triangles[it].vertices = std::array<VerIdx, 3>{iv1, iv2, iv3};
        m_triangulation.AddAdjacentTriangle(iv1, it);
        m_triangulation.AddAdjacentTriangle(iv2, it);
        m_triangulation.AddAdjacentTriangle(iv3, it);

        std::pair<TriIdx, TriIdx> ns;
        ns = m_triangulation.GetTriangles(IndexEdge(iv1, iv2));
        m_triangulation.ChangeNeighbor(ns.first, iv1, iv2, ns.second);
        m_triangulation.ChangeNeighbor(ns.second, iv1, iv2, ns.first);

        ns = m_triangulation.GetTriangles(IndexEdge(iv2, iv3));
        m_triangulation.ChangeNeighbor(ns.first, iv2, iv3, ns.second);
        m_triangulation.ChangeNeighbor(ns.second, iv2, iv3, ns.first);

        ns = m_triangulation.GetTriangles(IndexEdge(iv3, iv1));
        m_triangulation.ChangeNeighbor(ns.first, iv3, iv1, ns.second);
        m_triangulation.ChangeNeighbor(ns.second, iv3, iv1, ns.first);
        
        return it;
    }
    
    void RemoveOneTriangle(TriIdx it)
    {   
        for(VerIdx iv : m_triangulation.triangles[it].vertices)
            m_triangulation.RemoveAdjacentTriangle(iv, it);

        for(TriIdx neighbor : m_triangulation.triangles[it].neighbors){
            if(noNeighbor == neighbor) continue; 
            m_triangulation.ChangeNeighbor(neighbor, it, noNeighbor);
        }
        m_triangulation.triangles[it].Reset();
        m_removedTriangles.insert(it);
    }

    template <typename std::enable_if<Point::dim == 2, bool>::type = true>
    std::pair<TriIdx, TriIdx> TraversalTriangleAt(const Point & pos) const
    {
        auto out = std::make_pair(noNeighbor, noNeighbor);
        auto start = m_vertexIndexTree.NearestVertex(pos);
        auto triIdx = TraversalTriangles(start, pos);
        const auto & t  = m_triangulation.triangles[triIdx];
        const auto & p1 = m_triangulation.TrianglePoint(t, 0);
        const auto & p2 = m_triangulation.TrianglePoint(t, 1);
        const auto & p3 = m_triangulation.TrianglePoint(t, 2);
        auto loc = GetPointTriangleLocation(pos, p1, p2, p3);
        if(loc == PointTriangleLocation::Outside)
            throw std::runtime_error("No triangle was found at position");

        out.first = triIdx;
        if(loc < PointTriangleLocation::Inside && loc > PointTriangleLocation::Outside)
            out.second = m_triangulation.triangles[triIdx].neighbors[int(loc) - int(PointTriangleLocation::OnEdge1)];
        return out;
    }

    template <typename std::enable_if<Point::dim == 2, bool>::type = true>
    TriIdx TraversalTriangles(const VerIdx & start, const Point & pos) const
    {
        TriIdx curr = *(m_triangulation.vertices[start].triangles.begin());
        TriIdxFlatSet visited;
        bool found = false;
        while(!found){
            const auto & t = m_triangulation.triangles[curr];
            found = true;
            
            #ifdef GENERIC_UNIT_TEST
            auto offset = 0;
            #else
            auto offset = math::Random<index_t>(0, 2);
            #endif
            for(auto i = 0; i < 3; ++i){
                auto n = (i + offset) % 3;
                auto sVerPos = m_triangulation.TrianglePoint(t, n);
                auto eVerPos = m_triangulation.TrianglePoint(t, IndexTriangle::ccw(n));
                PointLineLocation edgeCheck = GetPointLineLocation(pos, sVerPos, eVerPos);
                if(edgeCheck == PointLineLocation::Right &&
                    t.neighbors[n] != noNeighbor &&
                    visited.insert(t.neighbors[n]).second){
                    found = false;
                    curr = t.neighbors[n];
                    break;
                }
            }
        }
        return curr;
    }

    void RebuildVertexIndexTree()
    {
        m_vertexIndexTree.Clear();
        for(VerIdx iv = 0; iv < m_triangulation.vertices.size(); ++iv)
            m_vertexIndexTree.AddOneVertex(m_triangulation.VerIdxPoint(iv), iv);
    }

    void ReallocateTriangulation()
    {
        VerIdxMap verIdxMap;
        TriIdxMap triIdxMap;
        ReallocateTriangulation(verIdxMap, triIdxMap);
    }

    bool ReallocateTriangulation(VerIdxMap & verIdxMap, TriIdxMap & triIdxMap)
    {
        if(!m_removedVertices.size() && !m_removedTriangles.size()) return false;
    
        verIdxMap.clear();    
        triIdxMap.clear();
        verIdxMap.insert(std::make_pair(noVertex, noVertex));
        triIdxMap.insert(std::make_pair(noNeighbor, noNeighbor));

        //points
        std::unordered_set<index_t> removedPoints;
        std::unordered_map<index_t, index_t> ptIdxMap;
        for(VerIdx iv : m_removedVertices)
            removedPoints.insert(m_triangulation.vertices[iv].index);
        
        for(index_t ip = 0, ipNew = 0; ip < m_triangulation.points.size(); ++ip){
            if(removedPoints.count(ip)) continue;
            ptIdxMap.insert(std::make_pair(ip, ipNew));
            m_triangulation.points[ipNew] = m_triangulation.points[ip];
            ipNew++;
        }
        m_triangulation.points.erase(m_triangulation.points.end() - removedPoints.size(), m_triangulation.points.end());

        //vertices
        for(VerIdx iv = 0, ivNew = 0; iv < m_triangulation.vertices.size(); ++iv){
            if(m_removedVertices.count(iv)) continue;
            verIdxMap.insert(std::make_pair(iv, ivNew));
            m_triangulation.vertices[ivNew] = m_triangulation.vertices[iv];
            ivNew++;
        }
        m_triangulation.vertices.erase(m_triangulation.vertices.end() - m_removedVertices.size(), m_triangulation.vertices.end());
        for(auto & vertex : m_triangulation.vertices){
            vertex.index = ptIdxMap.at(vertex.index);
        }
        
        //triangles
        for(TriIdx it = 0, itNew = 0; it < m_triangulation.triangles.size(); ++it){
            if(m_removedTriangles.count(it)) continue;
            triIdxMap.insert(std::make_pair(it, itNew));
            m_triangulation.triangles[itNew] = m_triangulation.triangles[it];
            itNew++;
        }
        m_triangulation.triangles.erase(m_triangulation.triangles.end() - m_removedTriangles.size(), m_triangulation.triangles.end());

        for(auto & vertex : m_triangulation.vertices){
            TriIdxSet triangles;
            for(auto it : vertex.triangles)
                triangles.insert(triIdxMap.at(it));
            std::swap(vertex.triangles, triangles);
        }

        for(auto & triangle : m_triangulation.triangles){
            std::array<VerIdx, 3> vertices;
            for(auto i = 0; i < 3; ++i)
                vertices[i] = verIdxMap.at(triangle.vertices[i]);
            std::swap(triangle.vertices, vertices);

            std::array<TriIdx, 3> neighbors;
            for(auto i = 0; i < 3; ++i)
                neighbors[i] = triIdxMap.at(triangle.neighbors[i]);
            std::swap(triangle.neighbors, neighbors);
        }

        //fixed edges
        EdgeSet fixedEdges;
        for(auto & e : m_triangulation.fixedEdges){
            fixedEdges.insert(IndexEdge(verIdxMap.at(e.v1()), verIdxMap.at(e.v2())));
        }
        std::swap(m_triangulation.fixedEdges, fixedEdges);

        m_removedTriangles.clear();
        m_removedVertices.clear();

        RebuildVertexIndexTree();
        return true;
    }

    /* Insert point into triangle: split into 3 triangles:
    *  - create 2 new triangles
    *  - re-use old triangle for the 3rd
    *                      v3
    *                    / | \
    *                   /  |  \ <-- original triangle (t)
    *                  /   |   \
    *          e3, n3 /    |    \ e2, n2
    *                /newT2|newT1\
    *               /      v      \
    *              /    __/ \__    \
    *             /  __/       \__  \
    *            / _/      t'     \_ \
    *          v1 ___________________ v2
    *                     e1, n1
    */
    VerIdx InsertPointInTriangle(const Point & pos, TriIdx it0, std::stack<TriIdx> & out)
    {    
        std::array<VerIdx, 3> vv = m_triangulation.triangles[it0].vertices;
        std::array<TriIdx, 3> nn = m_triangulation.triangles[it0].neighbors;

        VerIdx iv1 = vv[0], iv2 = vv[1], iv3 = vv[2];
        TriIdx in1 = nn[0], in2 = nn[1], in3 = nn[2];

        TriIdx it1 = AddOneTriangle();
        TriIdx it2 = AddOneTriangle();
        VerIdx iv  = AddOneVertex(pos, {it0, it1, it2});

        m_triangulation.triangles[it1].vertices = std::array<VerIdx, 3>{iv2, iv3, iv};
        m_triangulation.triangles[it1].neighbors = std::array<TriIdx, 3>{in2, it2, it0};

        m_triangulation.triangles[it2].vertices = std::array<VerIdx, 3>{iv3, iv1, iv};
        m_triangulation.triangles[it2].neighbors = std::array<TriIdx, 3>{in3, it0, it1};

        m_triangulation.triangles[it0].vertices = std::array<VerIdx, 3>{iv1, iv2, iv};
        m_triangulation.triangles[it0].neighbors = std::array<TriIdx, 3>{in1, it1, it2};

        m_triangulation.AddAdjacentTriangle(iv1, it2);
        m_triangulation.AddAdjacentTriangle(iv2, it1);
        m_triangulation.RemoveAdjacentTriangle(iv3, it0);
        m_triangulation.AddAdjacentTriangle(iv3, it1);
        m_triangulation.AddAdjacentTriangle(iv3, it2);

        m_triangulation.ChangeNeighbor(in2, it0, it1);
        m_triangulation.ChangeNeighbor(in3, it0, it2);

        out.push(it0);
        out.push(it1);
        out.push(it2);
        return iv;
    }

    /* Inserting a point on the edge between two triangles
    *    T1 (top)        v1
    *                   /|\
    *              n1 /  |  \ n4
    *               /    |    \
    *             /  T1' | Tnew1\
    *           v2-------v-------v4
    *             \ Tnew2| T2'  /
    *               \    |    /
    *              n2 \  |  / n3
    *                   \|/
    *   T2 (bottom)      v3
    */
    VerIdx InsertPointOnSharedEdge(const Point & pos, TriIdx it1, TriIdx it2, std::stack<TriIdx> & out)
    {
        auto & t1 = m_triangulation.triangles[it1];
        auto & t2 = m_triangulation.triangles[it2];
        index_t i  = t1.TnOpiVe(it2);
        VerIdx iv1 = t1.vertices[i];
        VerIdx iv2 = t1.vertices[IndexTriangle::ccw(i)];
        TriIdx in1 = t1.neighbors[i];
        TriIdx in4 = t1.neighbors[IndexTriangle::cw(i)];
        index_t j  = t2.TnOpiVe(it1);
        VerIdx iv3 = t2.vertices[j];
        VerIdx iv4 = t2.vertices[IndexTriangle::ccw(j)];
        TriIdx in3 = t2.neighbors[j];
        TriIdx in2 = t2.neighbors[IndexTriangle::cw(j)];

        TriIdx itNew1 = AddOneTriangle();
        TriIdx itNew2 = AddOneTriangle();
        VerIdx iv = AddOneVertex(pos, {it1, itNew2, it2, itNew1});

        m_triangulation.triangles[it1].vertices = std::array<VerIdx, 3>{iv1, iv2, iv};
        m_triangulation.triangles[it1].neighbors = std::array<TriIdx, 3>{in1, itNew2, itNew1};

        m_triangulation.triangles[it2].vertices = std::array<VerIdx, 3>{iv3, iv4, iv};
        m_triangulation.triangles[it2].neighbors = std::array<TriIdx, 3>{in3, itNew1, itNew2};

        m_triangulation.triangles[itNew1].vertices = std::array<VerIdx, 3>{iv1, iv, iv4};
        m_triangulation.triangles[itNew1].neighbors = std::array<TriIdx, 3>{it1, it2, in4};
        
        m_triangulation.triangles[itNew2].vertices = std::array<VerIdx, 3>{iv3, iv, iv2};
        m_triangulation.triangles[itNew2].neighbors = std::array<TriIdx, 3>{it2, it1, in2};

        m_triangulation.ChangeNeighbor(in4, it1, itNew1);
        m_triangulation.ChangeNeighbor(in2, it2, itNew2);
        m_triangulation.AddAdjacentTriangle(iv1, itNew1);
        m_triangulation.AddAdjacentTriangle(iv3, itNew2);
        m_triangulation.RemoveAdjacentTriangle(iv2, it2);
        m_triangulation.AddAdjacentTriangle(iv2, itNew2);
        m_triangulation.RemoveAdjacentTriangle(iv4, it1);
        m_triangulation.AddAdjacentTriangle(iv4, itNew1);

        out.push(it1);
        out.push(itNew2);
        out.push(it2);
        out.push(itNew1);
        return iv;
    }

    /* Inserting a point on the edge on boundary triangle
    *    T1 (top)        v1
    *                   /|\
    *              n1 /  |  \ n3
    *               /    |    \
    *             /  T'  | Tnew \
    *           v2-------v-------v3
    *                    n2
    */
    VerIdx InsertPointOnBoundaryEdge(const Point & pos, TriIdx it, const IndexEdge & e, std::stack<TriIdx> & out)
    {
        TriIdx itNew = AddOneTriangle();
        VerIdx iv = AddOneVertex(pos, {it, itNew});

        auto & t = m_triangulation.triangles[it];
        index_t ie = t.iEg(e);
        index_t i  = IndexTriangle::iEgOpiVe(ie);
        VerIdx iv1 = t.vertices[i];
        VerIdx iv2 = t.vertices[IndexTriangle::ccw(i)];
        VerIdx iv3 = t.vertices[IndexTriangle::cw(i)];
        TriIdx in1 = t.neighbors[i];
        TriIdx in2 = t.neighbors[IndexTriangle::ccw(i)];
        TriIdx in3 = t.neighbors[IndexTriangle::cw(i)];
        assert(in2 == noNeighbor);

        m_triangulation.triangles[it].vertices = std::array<VerIdx, 3>{iv1, iv2, iv};
        m_triangulation.triangles[it].neighbors = std::array<TriIdx, 3>{in1, in2, itNew};
        m_triangulation.triangles[itNew].vertices = std::array<VerIdx, 3>{iv1, iv, iv3};
        m_triangulation.triangles[itNew].neighbors = std::array<TriIdx, 3>{it, in2, in3};

        m_triangulation.ChangeNeighbor(in3, it, itNew);
        m_triangulation.AddAdjacentTriangle(iv1, itNew);
        m_triangulation.RemoveAdjacentTriangle(iv3, it);
        m_triangulation.AddAdjacentTriangle(iv3, itNew);

        out.push(it);
        out.push(itNew);
        return iv;
    }

    void Delaunay(std::stack<TriIdx> & triangles, VerIdx iv)
    {
        auto fixed = [this](TriIdx it, VerIdx iv)
        {
            auto e = m_triangulation.triangles[it].VeOpEg(iv);
            return m_triangulation.fixedEdges.count(e);
        };

        while(!triangles.empty()){
            TriIdx it = triangles.top();
            triangles.pop();

            const auto & t = m_triangulation.triangles[it];
            TriIdx itOp = t.VeOpTn(iv);
            if(itOp == noNeighbor) continue;
            if(fixed(it, iv)) continue;
            if(!Utility::isLocallyDelaunay(m_triangulation, it, itOp)){
                m_triangulation.FlipEdge(it, itOp);
                triangles.push(it);
                triangles.push(itOp);
            }
        }
    }

    bool TryDecayOneEdge(const IndexEdge & e)
    {
        VerIdx iv1  = e.v1();
        VerIdx iv2  = e.v2();
        const auto & p1 = m_triangulation.VerIdxPoint(iv1);
        const auto & p2 = m_triangulation.VerIdxPoint(iv2);
        Point p = (p1 + p2) * 0.5;
        auto modifiable = [&, this](const TriIdx it, const Point & p)
        {
            if(noNeighbor == it) return true;
            const auto & t = m_triangulation.triangles[it];

            VerIdx ivOp = t.EgOpVe(e);
            TriIdx itn1 = t.Tn(ivOp, iv1);
            TriIdx itn2 = t.Tn(ivOp, iv2);
            if(noNeighbor != itn1){
                const auto & tn1 = m_triangulation.triangles[itn1];
                auto [iu, iw] = tn1.VeOpVes(iv1);
                const auto & u = m_triangulation.VerIdxPoint(iu);
                const auto & w = m_triangulation.VerIdxPoint(iw);
                if(!Triangle::isCCW(p, w, u)) return false;                
            }

            if(noNeighbor != itn2){
                const auto & tn2 = m_triangulation.triangles[itn2];
                auto [iu, iw] = tn2.VeOpVes(iv2);
                const auto & u = m_triangulation.VerIdxPoint(iu);
                const auto & w = m_triangulation.VerIdxPoint(iw);
                if(!Triangle::isCCW(p, w, u)) return false;
            }
            return true;
        };
        auto [it1, it2] = m_triangulation.GetTriangles(e);
        if(!modifiable(it1, p)) return false;
        if(!modifiable(it2, p)) return false;

        DecayToVertex(e);
        return true;
    }

    void DecayToVertex(const IndexEdge & e)
    {
        VerIdx iv1  = e.v1();
        VerIdx iv2  = e.v2();
        const auto & p1 = m_triangulation.VerIdxPoint(iv1);
        const auto & p2 = m_triangulation.VerIdxPoint(iv2);
        Point p = (p1 + p2) * 0.5;
        m_triangulation.points[m_triangulation.vertices[iv1].index] = p;
        auto removeOneTriangle = [&, this](TriIdx it) mutable
        {
            if(noNeighbor == it) return;
            const auto & t = m_triangulation.triangles[it];

            VerIdx ivOp = t.EgOpVe(e);
            TriIdx itn1 = t.Tn(ivOp, iv1);
            TriIdx itn2 = t.Tn(ivOp, iv2);

            RemoveOneTriangle(it);
            m_triangulation.ChangeNeighbor(itn1, iv1, ivOp, itn2);
            m_triangulation.ChangeNeighbor(itn2, iv2, ivOp, itn1);
            if(noNeighbor != itn1 && !Utility::isCCW(m_triangulation, itn1))
                m_triangulation.triangles[itn1].Reverse();
            if(noNeighbor != itn2 && !Utility::isCCW(m_triangulation, itn2))
                m_triangulation.triangles[itn2].Reverse();
        };
        auto [it1, it2] = m_triangulation.GetTriangles(e);
        removeOneTriangle(it1);
        removeOneTriangle(it2);
        for(auto it : m_triangulation.vertices[iv2].triangles){
            auto & t = m_triangulation.triangles[it];
            t.vertices[t.iVe(iv2)] = iv1;
            auto & v = m_triangulation.vertices[iv1];
            v.triangles.insert(it);
        }

        RemoveOneVertex(iv2);

        std::list<IndexEdge> toRemove, toInsert;
        for(const auto & edge : m_triangulation.fixedEdges){
            if(edge.hasVertex(iv2)){
                toRemove.push_back(edge);
                if(edge != e){
                    if(edge.v1() == iv2) toInsert.emplace_back(IndexEdge{iv1, edge.v2()});
                    else toInsert.emplace_back(IndexEdge{edge.v1(), iv1});
                }
            }
        }

        for(const auto & edge : toRemove)
            m_triangulation.fixedEdges.erase(edge);
        for(const auto & edge : toInsert)
            m_triangulation.fixedEdges.insert(edge);
    }

    void Clear()
    {
        m_triangulation.Clear();
        m_vertexIndexTree.Clear();
        m_removedTriangles.clear();
        m_removedVertices.clear();
    }
};

}//namespace tri
}//namespace geometry
}//namespace generic