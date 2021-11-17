#ifndef GENERIC_GEOMETRY_TRI_TRIANGULATIONREFINEMENT_HPP
#define GENERIC_GEOMETRY_TRI_TRIANGULATIONREFINEMENT_HPP
#include "TriangleEvaluator.hpp"
#include "Triangulator.hpp"
#include "Utility.hpp"
#include <unordered_map>
#include <unordered_set>
#include <optional>
#include <tuple>
#include <queue>
#include <stack>
namespace generic  {
namespace geometry {
namespace tri      {

template <typename num_type>
class DelaunayRefinement2D
{
protected:
    using float_t = float_type<num_type>;
    using Point = Point2D<num_type>;
    using PointF = Point2D<float_t>;
    using Utility = TriangulationUtility<Point2D<num_type> >;
    using EdgeSet = typename Triangulation<Point>::EdgeSet;

    struct Mark {
        bool fe = false;
        bool ft = false;
        bool onEdge = false;
        VerIdx ive = noVertex;
        TriIdx its = noNeighbor;
        IndexEdge e = IndexEdge(noVertex, noVertex);
        IndexEdge se = IndexEdge(noVertex, noVertex);//test, bwu
    };

    struct State {
        Mark prev;
        Mark curr;
    };

    struct Paras {
        float_t minAlpha = 0;
        float_t maxBoundSq = std::numeric_limits<float_t>::max();
        float_t maxOcBoundSq = std::numeric_limits<float_t>::max();
        float_t minEdgeLenSq = 0;
        float_t maxEdgeLenSq = std::numeric_limits<float_t>::max();
    };

    State state;
    Paras paras;
    Triangulation<Point> & tri;
    TriangulationOperator<Point> op;
public:
    explicit DelaunayRefinement2D(Triangulation<Point> & t) : tri(t), op(t){}
    virtual ~DelaunayRefinement2D() = default;

    const TriangulationOperator<Point> & GetOp() const { return op; }//test, bwu
    void SetParas(float_t alpha, num_type minEdgeLen, num_type maxEdgeLen)
    {
        auto limitAlpha = LimitAlpha();
        if(alpha < 0) alpha = limitAlpha;
        paras.minAlpha = math::LE(alpha, limitAlpha) ? alpha : limitAlpha;
        float_t bound = 0.5 / std::sin(paras.minAlpha);
        float_t ocBound = 0.5 / std::sin(0.5 * paras.minAlpha);
        paras.maxBoundSq = bound * bound;
        paras.maxOcBoundSq = ocBound * ocBound;
        paras.minEdgeLenSq = minEdgeLen * minEdgeLen;
        paras.maxEdgeLenSq = maxEdgeLen * maxEdgeLen;
        UpdateState();
    }
    virtual void MergeShortestEdge(){}//test, bwu
    virtual void UpdateState() = 0;
    virtual void Refine(index_t step) = 0;
    virtual float_t LimitAlpha() const = 0;
    virtual std::optional<Point2D<num_type> > CurrVertexPoint() const
    {
        if(state.curr.ft){
            if(isBadlySkinnyTriangle(state.curr.its))
                return GetOffCenterPoint(state.curr.its).template Cast<num_type>();
            else return GetCircumCenterPoint(state.curr.its).template Cast<num_type>();
        }
        return std::nullopt;
    }

    virtual std::optional<Segment2D<num_type> > CurrEncroachedEdge() const
    {
        if(state.curr.fe) return Utility::GetSegment(tri, state.curr.e);
        return std::nullopt;
    }

    virtual std::optional<Triangle2D<num_type> > CurrSkinnyTriangle() const
    {
        if(state.curr.ft) return Utility::GetTriangle(tri, state.curr.its);
        return std::nullopt;
    }

    virtual std::optional<Segment2D<num_type> > CurrShortestEdge() const//test, bwu
    {
        if(stillExist(state.curr.se))
            return Utility::GetSegment(tri, state.curr.se);
        return std::nullopt;
    }

protected:
    PointF GetCircumCenterPoint(const TriIdx & it) const
    {
        return Utility::GetCircumCenter(tri, it);
    }

    PointF GetOffCenterPoint(const TriIdx & it) const
    {
        PointF cc = GetCircumCenterPoint(it);
        IndexEdge e = Utility::GetMinLenEdge(tri, it);
        PointF p1 = Utility::GetVertexPoint(tri, e.v1()).template Cast<float_t>();
        PointF p2 = Utility::GetVertexPoint(tri, e.v2()).template Cast<float_t>();
        auto dir = (p2 - p1);
        auto mid = (p1 + p2) * 0.5;
        auto norm = geometry::Normalize(PointF(-dir[1], dir[0]));
        auto alpha = std::max(0.5 * paras.minAlpha, TriangleEvaluator<Point>::MinimumAngle(tri, it));
        auto t = 0.5 * geometry::Distance(p1, p2) / std::tan(alpha);

        PointF oc = mid + norm * t;
        if(geometry::DotProduct(cc - p1, oc - p1) > 0) return oc;
        else return mid - norm * t;
    }

    VerIdx InsertOnePointInTriangulation(const Point & p)
    {
        VerIdx iv;
        std::stack<TriIdx> tris;
        auto [it1, it2] = op.TraversalTriangleAt(p);
        if(noNeighbor == it2){
            auto loc = Utility::GetPointTriangleLocation(tri, it1, p);
            GENERIC_ASSERT(loc != PointTriangleLocation::Outside)
            if(PointTriangleLocation::Inside == loc)
                iv = op.InsertPointInTriangle(p, it1, tris);
            else {
                auto edge = tri.triangles[it1].Edge(static_cast<int>(loc) - 1);
                iv = op.InsertPointOnBoundaryEdge(p, it1, edge, tris);
            }
        }
        else iv = op.InsertPointOnSharedEdge(p, it1, it2, tris);
        op.Delaunay(tris, iv);
        return iv;
    }

    void RemoveOneVertexFromTriangulation(const VerIdx & iv, std::vector<VerIdx> & free)
    {
        Utility::GetNeighborVerticesCCW(tri, iv, free);

        auto its = tri.vertices[iv].triangles;
        for(TriIdx it : its)
            op.RemoveOneTriangle(it);
        op.RemoveOneVertex(iv);
    }

    void RemoveOneVertexFromTriangulation(const VerIdx & iv, EdgeSet & free)
    {
        free.clear();
        IndexVertex & vertex = tri.vertices[iv];
        auto its = vertex.triangles;
        for(TriIdx it : its)
            free.insert(tri.triangles[it].VeOpEg(iv));
        for(TriIdx it : its)
            op.RemoveOneTriangle(it);
        op.RemoveOneVertex(iv);
    }

    void ReTriangulate(const std::vector<VerIdx> & free)
    {
        std::list<TriIdx> triangles;
        ReTriangulate(free, triangles);
    }

    void ReTriangulate(const std::vector<VerIdx> & free, std::list<TriIdx> & triangles)
    {
        Triangulation<Point> localTriangulation;
        ConvexTriangulator2D<num_type> localTriangulator(localTriangulation);
        auto xGetter = [this](const VerIdx & iv){ return tri.points[tri.vertices[iv].index][0]; };
        auto yGetter = [this](const VerIdx & iv){ return tri.points[tri.vertices[iv].index][1]; };
        localTriangulator.Triangulate(free.begin(), free.end(), xGetter, yGetter);

        for(const auto & t : localTriangulation.triangles){
            auto it = op.AddOneTriangle(free[t.vertices[0]], free[t.vertices[1]], free[t.vertices[2]]);
            triangles.push_back(it);
        }
    }

    void ReTriangulate(const EdgeSet & free)
    {
        std::list<TriIdx> triangles;
        ReTriangulate(free, triangles);
    }


    void ReTriangulate(const EdgeSet & free, std::list<TriIdx> & triangles)
    {
        if(0 == free.size()) return;

        std::vector<Point> points;
        std::unordered_map<VerIdx, index_t> viMap;
        std::unordered_map<index_t, VerIdx> ivMap;
        std::list<std::pair<size_t, size_t> > edges;
        auto addVertex = [&](const VerIdx & iv) mutable
        {
            if(viMap.count(iv)) return;
            index_t i = points.size();
            ivMap.insert(std::make_pair(i, iv));
            viMap.insert(std::make_pair(iv, i));
            points.push_back(Utility::GetVertexPoint(tri, iv));
        };

        for(const auto & e : free){
            addVertex(e.v1());
            addVertex(e.v2());
            edges.push_back(std::make_pair(viMap.at(e.v1()), viMap.at(e.v2())));
        }

        Triangulation<Point> localTriangulation;
        Triangulator2D<num_type> localTriangulator(localTriangulation);
        localTriangulator.InsertVertices(points.begin(), points.end(), [](const Point & p){ return p[0]; }, [](const Point & p){ return p[1]; });
        localTriangulator.InsertEdges(edges.begin(), edges.end(), [](const std::pair<size_t, size_t> & e){ return e.first; }, [](const std::pair<size_t, size_t> & e){ return e.second; });
        localTriangulator.EraseOuterTriangles();

        for(const auto & t : localTriangulation.triangles){
            auto it = op.AddOneTriangle(ivMap[t.vertices[0]], ivMap[t.vertices[1]], ivMap[t.vertices[2]]);
            triangles.push_back(it);
        }
    }
    
    bool stillExist(const IndexEdge & e) const
    {
        if(noVertex == e.v1()) return false;
        if(noVertex == e.v2()) return false;
        const auto & v1 = tri.vertices[e.v1()];
        const auto & v2 = tri.vertices[e.v2()];
        return IndexVertex::isShareEdge(v1, v2);
    }

    Point SplitEdgeAtCircularShell(const IndexEdge & e)//for float_t,//test, bwu
    {
        const auto & p1 = tri.VerIdxPoint(e.v1());
        const auto & p2 = tri.VerIdxPoint(e.v2());

        size_t i = 0;
        float_t low, high, mid = geometry::Distance(p1, p2) * 0.5;
        for(; i < sizeof(num_type) * 8; ++i){
            low  = std::pow(2, i);
            high = std::pow(2, i + 1);
            if(math::GT(mid, low) && math::LE(mid, high)) break;
        }
        float_t t = std::pow(2, i + 1) / (2.0 * mid);
        return p1 + (p2 - p1) * t;
    }

    Point SplitEdgeAtMidPoint(const IndexEdge & e)
    {
        const auto & p1 = tri.VerIdxPoint(e.v1());
        const auto & p2 = tri.VerIdxPoint(e.v2());
        auto mid = (p1 + p2) * 0.5;
        return mid;
    }

    void FindEncroachedCircumCenters(const VerIdxSet & circumcenters, const IndexEdge & e, std::list<VerIdx> & encroached)
    {
        TriIdxSet visited;
        auto [it1, it2] = tri.GetTriangles(e);
        auto circle = Utility::GetDiametralCircle(tri, e);
        TraversalEncroachedCircumCentersAt(circumcenters, e, it1, circle, visited, encroached);
        TraversalEncroachedCircumCentersAt(circumcenters, e, it2, circle, visited, encroached);
    }

    void TraversalEncroachedCircumCentersAt(const VerIdxSet & circumcenters, const IndexEdge & e, const TriIdx it, const Circle<float_t> & c, TriIdxSet & visited, std::list<VerIdx> & encroached)
    {
        if(noNeighbor == it) return;
        if(visited.count(it)) return;

        visited.insert(it);
        const auto & t = tri.triangles[it];
        VerIdx iv = t.EgOpVe(e);
        
        bool contains = c.Contains(tri.VerIdxPoint(iv));
        if(contains && circumcenters.count(iv))
            encroached.push_back(iv);

        auto iu = t.NextVertex(iv, WindingDirection::Clockwise);
        auto iw = t.NextVertex(iv, WindingDirection::CounterClockwise);
        if(c.Contains(tri.VerIdxPoint(iu))) TraversalEncroachedCircumCentersAt(circumcenters, IndexEdge(iv, iu), t.Tn(iv, iu), c, visited, encroached);
        if(c.Contains(tri.VerIdxPoint(iw))) TraversalEncroachedCircumCentersAt(circumcenters, IndexEdge(iv, iw), t.Tn(iv, iw), c, visited, encroached);
    }

    bool FindFirstSkinnyTriangle(TriIdx & it) const
    {
        index_t size = tri.triangles.size();
        index_t start = math::Random(index_t(0), size - 1);
        it = start;
        float_t ratioSq;
        do{
            if(!op.isTriangleRemoved(it)){
                if(isSkinnyTriangle(it, ratioSq))
                    return true;
            }
            it = (it + 1) % size;
        } while(it != start);
        return false;
    }

    bool isSkinnyTriangle(const TriIdx & it) const
    {
        float_t ratioSq;
        return isSkinnyTriangle(it, ratioSq);
    }

    bool isSkinnyTriangle(const TriIdx & it, float_t & ratioSq) const
    {
        // float_t aveEdgeLen = Utility::GetAverageEdgeLength(tri, it);//test, bwu
        // if(aveEdgeLen * aveEdgeLen < paras.minEdgeLenSq) return false;//test, bwu
        ratioSq = Utility::CircumRadius2ShortestEdgeRatioSq(tri, it);
        return ratioSq > paras.maxBoundSq;
    }

    bool isBadlySkinnyTriangle(const TriIdx & it) const
    {
        float_t ratioSq = Utility::CircumRadius2ShortestEdgeRatioSq(tri, it);
        return ratioSq > paras.maxOcBoundSq;
    }

    bool isOverSizedTriangle(const TriIdx & it) const
    {
        auto lenSq = Utility::GetMaxEdgeLenSq(tri, it);
        return lenSq > paras.maxEdgeLenSq;
    }

    bool isTinySizedTriangle(const TriIdx & it) const
    {
        auto lenSq = Utility::GetMaxEdgeLenSq(tri, it);
        return lenSq < paras.minEdgeLenSq;
    }
};

template <typename num_type>
class JonathanRefinement2D : public DelaunayRefinement2D<num_type>
{
    using VerIdxMap = std::unordered_map<VerIdx, VerIdx>;
    using TriIdxMap = std::unordered_map<TriIdx, TriIdx>;
    using Base = DelaunayRefinement2D<num_type>;
    using float_t = typename Base::float_t;
    using EdgeSet = typename Base::EdgeSet;
    using PointF = typename Base::PointF;
    using Point = typename Base::Point;
    using Base::state;
    using Base::paras;
    using Base::tri;
    using Base::op;

    struct EncroachedEdge
    {
        IndexEdge edge;
        float_t weight;
    };

    struct EncroachedEdgeCmp
    {
        bool operator() (const EncroachedEdge & e1, const EncroachedEdge & e2) const noexcept
        {
            return e1.weight < e2.weight;
        }
    };
    using EncroachedEdgeQueue = std::priority_queue<EncroachedEdge, std::vector<EncroachedEdge>, EncroachedEdgeCmp>;
    
    struct SkinnyTriangle
    {
        TriIdx index;
        float_t weight;
    };

    struct SkinnyTriangleCmp
    {
        bool operator() (const SkinnyTriangle & t1, const SkinnyTriangle & t2) const noexcept
        {
            return t1.weight < t2.weight;
        }
    };

    using SkinnyTriangleQueue = std::priority_queue<SkinnyTriangle, std::vector<SkinnyTriangle>, SkinnyTriangleCmp>;

    class DiametralLens
    {
    public:
        DiametralLens(Point start, Point end) : s(start), e(end)
        {
            auto theta = math::Rad(30);
            auto cosT = std::cos(theta);
            auto sinT = std::sin(theta);
            auto sf = s.template Cast<float_t>();
            auto ef = e.template Cast<float_t>();
            auto vf = (ef - sf) / std::sqrt(3);
            auto cw  = sf + PointF(vf[0] * cosT + vf[1] * sinT, vf[1] * cosT - vf[0] * sinT);
            auto ccw = sf + PointF(vf[0] * cosT - vf[1] * sinT, vf[1] * cosT + vf[0] * sinT);
            cl = geometry::CircumCircle( cw, sf, ef);
            cr = geometry::CircumCircle(ccw, sf, ef); 
        }

        bool Contains(const Point & p) const
        {
            using segment_t = segment_type<Point>;
            auto loc = geometry::GetPointLineLocation(p, s, e);
            if(geometry::PointLineLocation::Left == loc) return cr.Contains(p);
            else if(geometry::PointLineLocation::Right == loc) return cl.Contains(p);
            else if(geometry::Contains(segment_t(s, e), p, true))return true;
            return false;
        }
    private:
        Point s, e;//start,end
        Circle<float_t> cl, cr;//left,right
    };

    struct EdgeSplitResult
    {
        VerIdx split;
        IndexEdge e1;
        IndexEdge e2;
    };

    std::queue<IndexEdge> m_queueSE;//shortest edge, for test, bwu
public:
    explicit JonathanRefinement2D(Triangulation<Point> & t)
     : DelaunayRefinement2D<num_type>(t)
    {
    }

    void UpdateShortestEdge()//test, bwu
    {
        for(auto it = 0; it < tri.triangles.size(); ++it){
            if(op.isTriangleRemoved(it)) continue;
            const auto & t = tri.triangles[it];
            for(auto ie = 0; ie < 3; ++ie){
                auto e = t.Edge(ie);
                auto distSq = Base::Utility::GetEdgeLenSq(tri, e);
                if(distSq < paras.minEdgeLenSq)
                    m_queueSE.push(std::move(e));
            }
        }
        if(!m_queueSE.empty())
            state.curr.se = m_queueSE.front();
    }

    float_t LimitAlpha() const { return math::Rad(30); }

    void MergeShortestEdge()
    {
        while(!m_queueSE.empty()){
            auto e = m_queueSE.front();
            m_queueSE.pop();
            if(!Base::stillExist(e)) continue;
            auto distSq = Base::Utility::GetEdgeLenSq(tri, e);
            if(distSq > paras.minEdgeLenSq) continue;
            op.TryDecayOneEdge(e);
        }
        UpdateShortestEdge();
    }

    void UpdateState() {
        m_queueE = EncroachedEdgeQueue{};
        m_queueT = SkinnyTriangleQueue{};
        EnqueueEncroachedEdges();
    }

    void SplitEncroachedEdges(index_t & maxStep)
    {
        while(maxStep && DequeueInvalidEdges()){
            auto e = m_queueE.top();
            m_queueE.pop();

            auto res = SplitEncroachedEdge(e.edge);
            std::list<TriIdx> triangles;
            RemoveEncroachedCircumCenters(e.edge, triangles);
            for(auto it : triangles){
                if(tri.vertices[res.split].triangles.count(it)) continue;
            }
            NewVertex(res.split, [](const TriIdx it){ return false; });
            if(isEncroached(res.e1)) m_queueE.push({res.e1, WE(res.e1)});
            if(isEncroached(res.e2)) m_queueE.push({res.e2, WE(res.e2)});
            maxStep--;
        }
        if(0 == m_queueE.size())
            EnqueueSkinnyTriangles();
    }

    void RemoveSkinnyTriangles(index_t & maxStep)
    {
        auto func = [this](const TriIdx it) { return Base::isOverSizedTriangle(it); };
        while(maxStep && DequeueInvalidTriangles()){
            auto t = m_queueT.top();
            m_queueT.pop();

            std::list<IndexEdge> edges;
            PointF cc = Base::Utility::GetCircumCenter(tri, t.index);
            GetEncroachedEdges(cc.template Cast<num_type>(), edges);
            if(edges.empty()){
                VerIdx iv = InsertCircumCenter(cc);
                NewVertex(iv, func);
            }
            else {
                auto minSq = Base::Utility::GetMinEdgeLenSq(tri, t.index);
                for(const auto & e : edges){
                    if(func(t.index) || SplitPermitted(e, minSq))
                        m_queueE.push({e, WE(e)});
                }
                if(!m_queueE.empty()){
                    m_queueT.push({t.index, WT(t.index)});
                    SplitEncroachedEdges(paras.minAlpha, func);
                }
            }
            maxStep--;
        }
    }

    void Refine(index_t maxStep)
    {
        SplitEncroachedEdges(maxStep);
        RemoveSkinnyTriangles(maxStep);
    }

    void ReallocateTriangulation()
    {
        VerIdxMap verIdxMap;
        TriIdxMap triIdxMap;
        if(op.ReallocateTriangulation(verIdxMap, triIdxMap)){
            //remap circumcenters
            VerIdxSet circumcenters;
            for(auto & iv : m_circumcenters)
                circumcenters.insert(verIdxMap.at(iv));
            std::swap(m_circumcenters, circumcenters);
        }
    }

    std::optional<Segment2D<num_type> > CurrEncroachedEdge() const
    {
        DequeueInvalidEdges();
        if(!m_queueE.empty())
            return Base::Utility::GetSegment(tri, m_queueE.top().edge);
        return std::nullopt;
    }

    std::optional<Triangle2D<num_type> > CurrSkinnyTriangle() const
    {
        DequeueInvalidTriangles();
        if(!m_queueT.empty())
            return Base::Utility::GetTriangle(tri, m_queueT.top().index);
        return std::nullopt;
    }

private:
    void EnqueueEncroachedEdges()
    {
        for(const auto & e : tri.fixedEdges){
            if(isEncroached(e))
                m_queueE.push({e, WE(e)});
        }
    }

    bool DequeueInvalidEdges() const
    {
        while(!m_queueE.empty() &&
              !Base::stillExist(m_queueE.top().edge)){
            m_queueE.pop();
        }
        return !m_queueE.empty();
    }

    void EnqueueSkinnyTriangles()
    {
        for(TriIdx it = 0; it < tri.triangles.size(); ++it){
            if(op.isTriangleRemoved(it)) continue;
            if(Base::isSkinnyTriangle(it)) m_queueT.push({it, WT(it)});
            else if(Base::isOverSizedTriangle(it)) m_queueT.push({it, WT(it)});
        }
    }

    bool DequeueInvalidTriangles() const
    {
        while(!m_queueT.empty() &&
              (op.isTriangleRemoved(m_queueT.top().index) ||
              (Base::isTinySizedTriangle(m_queueT.top().index)) ||
              (!Base::isSkinnyTriangle(m_queueT.top().index) &&
               !Base::isOverSizedTriangle(m_queueT.top().index)))){
            m_queueT.pop();
        }
        return !m_queueT.empty();
    }

    bool isTinyEdge(const IndexEdge & e) const
    {
        return Base::Utility::GetEdgeLenSq(e) < paras.minEdgeLenSq;
    }

    bool isEncroached(const IndexEdge & e) const
    {
        DiametralLens dl(tri.VerIdxPoint(e.v1()), tri.VerIdxPoint(e.v2()));
        auto [it1, it2] = tri.GetTriangles(e);

        auto ehcroached = [this, &e, &dl](TriIdx it)
        {
            const auto & t = tri.triangles[it];
            VerIdx iv = t.EgOpVe(e);
            return dl.Contains(tri.VerIdxPoint(iv));
        };

        if(noNeighbor != it1 && ehcroached(it1)) return true;
        if(noNeighbor != it2 && ehcroached(it2)) return true;
        return false;
    }

    bool isEncroached(const IndexEdge & e, const VerIdx iv) const
    {
        DiametralLens dl(tri.VerIdxPoint(e.v1()), tri.VerIdxPoint(e.v2()));
        return dl.Contains(tri.VerIdxPoint(iv));
    }

    bool SplitPermitted(const IndexEdge & e, const num_type lenSq) const//todo, float_type handle
    {
        float_t tolerance = 0.01;
        float_t p = 0.5 * std::log(lenSq) / std::log(2);
        if(lenSq > 4 && !math::EQ(std::round(p), p, tolerance)) return true;

        auto cluster1 = tri.ConcentricFixedEdges(e.v1());
        auto cluster2 = tri.ConcentricFixedEdges(e.v2());
        GENERIC_ASSERT(cluster1.size() > 0)
        GENERIC_ASSERT(cluster2.size() > 0)
        if(1 == cluster1.size() && 1 == cluster2.size()) return true;
        if(2 <= cluster1.size() && 2 <= cluster2.size()) return true;

        const auto & cluster = (2 <= cluster1.size()) ? cluster1 : cluster2;
        num_type minSq = std::numeric_limits<num_type>::max();
        for(const auto & edge : cluster){
            if(edge == e) continue;
            auto tmp = Base::Utility::GetEdgeLenSq(tri, edge);
            if(tmp < minSq) minSq = tmp;
        }
        if(math::LT(minSq, lenSq)) return true;
        if(math::GE(float_t(minSq / 4.0), float_t(lenSq))) return true;
        return false;
    }

    void GetEncroachedEdges(const Point & p, std::list<IndexEdge> & edges)
    {
        //todo enhance to avoid iteration
        edges.clear();
        for(const auto & e : tri.fixedEdges){
            DiametralLens dl(tri.VerIdxPoint(e.v1()), tri.VerIdxPoint(e.v2()));
            if(dl.Contains(p)) edges.push_back(e);
        }
    }

    template <typename Func>
    void SplitEncroachedEdges(float_t theta, Func && func)
    {
        while(DequeueInvalidEdges()){
            auto e = m_queueE.top();
            m_queueE.pop();
            if(Base::Utility::GetEdgeLenSq(tri, e.edge) < paras.minEdgeLenSq){
                continue;//test, bwu
            }

            auto res = SplitEncroachedEdge(e.edge);
            std::list<TriIdx> triangles;
            RemoveEncroachedCircumCenters(e.edge, triangles);
            for(auto it : triangles){
                if(tri.vertices[res.split].triangles.count(it)) continue;
                if(func(it)) m_queueT.push({it, WT(it)});
                else if(Base::Utility::GetMinimumAngle(tri, it) < theta) m_queueT.push({it, WT(it)});
            }
            NewVertex(res.split, func);
            if(isEncroached(res.e1)) m_queueE.push({res.e1, WE(res.e1)});
            if(isEncroached(res.e2)) m_queueE.push({res.e2, WE(res.e2)});
        }
    }

    template <typename Func>
    void NewVertex(VerIdx iv, Func && func)
    {
        const auto & vertex = tri.vertices[iv];
        for(auto it : vertex.triangles){
            const auto & triangle = tri.triangles[it];
            auto e = triangle.VeOpEg(iv);
            if(tri.hasFixedEdge(e) && isEncroached(e, iv)) m_queueE.push({e, WE(e)});
            else if(func(it)) m_queueT.push({it, WT(it)});
            else if(Base::Utility::GetMinimumAngle(tri, it) < paras.minAlpha) m_queueT.push({it, WT(it)});
        }
    }

    EdgeSplitResult SplitEncroachedEdge(const IndexEdge & e)
    {
        auto split = GetEdgeSplitPosition(e);//test, bwu

        VerIdx iv;
        std::stack<TriIdx> tris;
        auto [it, itOp] = tri.GetTriangles(e);
        if(itOp == noNeighbor)
            iv = op.InsertPointOnBoundaryEdge(split, it, e, tris);
        else iv = op.InsertPointOnSharedEdge(split, it, itOp, tris);
        op.Delaunay(tris, iv);

        EdgeSplitResult res{ iv, IndexEdge(e.v1(), iv), IndexEdge(iv, e.v2()) };

        if(tri.fixedEdges.count(e)){
            tri.fixedEdges.erase(e);
            tri.fixedEdges.insert(res.e1);
            tri.fixedEdges.insert(res.e2);
        }
        
        return res;
    }

    Point GetEdgeSplitPosition(const IndexEdge & e)
    {
        std::list<IndexEdge> edges;
        const auto & p1 = tri.VerIdxPoint(e.v1());
        const auto & p2 = tri.VerIdxPoint(e.v2());
        edges = tri.ConcentricFixedEdges(e.v1());
        if(edges.size() > 1){
            for(const auto & edge : edges){
                if(edge == e) continue;
                const auto & p = edge.v1() == e.v1() ? tri.VerIdxPoint(edge.v2()) : tri.VerIdxPoint(edge.v1());
                if(Angle(p, p1, p2) < math::pi_half)
                    return Base::SplitEdgeAtCircularShell(e);
            }
        }

        edges = tri.ConcentricFixedEdges(e.v2());
        if(edges.size() > 1){
            for(const auto & edge : edges){
                if(edge == e) continue;
                const auto & p = edge.v1() == e.v2() ? tri.VerIdxPoint(edge.v2()) : tri.VerIdxPoint(edge.v1());
                if(Angle(p, p2, p1) < math::pi_half)
                    return Base::SplitEdgeAtCircularShell(e);
            }
        }

        return Base::SplitEdgeAtMidPoint(e);
    }

    void RemoveEncroachedCircumCenters(const IndexEdge & e, std::list<TriIdx> & triangles)
    {
        EdgeSet free;//test, bwu
        // std::vector<VerIdx> free;//some problem here
        std::list<VerIdx> encroached;
        Base::FindEncroachedCircumCenters(m_circumcenters, e, encroached);
        for(auto iv : encroached){
            Base::RemoveOneVertexFromTriangulation(iv, free);
            Base::ReTriangulate(free, triangles);
            m_circumcenters.erase(iv);
        }
    }

    VerIdx InsertCircumCenter(const PointF & cc)
    {
        VerIdx iv = Base::InsertOnePointInTriangulation(cc.template Cast<typename Point::coor_t>());
        m_circumcenters.insert(iv);
        return iv;
    }

    float_t WE(const IndexEdge & e) const
    {
        return Base::Utility::GetEdgeLenSq(tri, e);
    }

    float_t WT(TriIdx it) const
    {
        // return Base::Utility::CircumRadius2ShortestEdgeRatioSq(tri, it);//test, bwu
        return Base::Utility::GetTriangleArea(tri, it);
    }

    bool isCapTriangle(TriIdx it) const
    {
        float_t threshold = 1e3;
        auto angles = Base::Utility::GetInteriorAngles(tri, it);
        std::sort(angles.begin(), angles.end(), std::less<float_t>());
        return (angles[2] / angles[1] > threshold) && (angles[2] / angles[0] > threshold);
    }

    void RemoveCapTriangle(TriIdx it)
    {
        auto e = Base::Utility::GetMaxLenEdge(tri, it);
        auto iv1 = tri.triangles[it].EgOpVe(e);
        auto iv2 = SplitEdge(e);
        op.TryDecayOneEdge(IndexEdge(iv1, iv2));
    }

    VerIdx SplitEdge(const IndexEdge & e)
    {
        std::stack<TriIdx> outs;
        auto [it1, it2] = tri.GetTriangles(e);
        Point pos = Base::SplitEdgeAtMidPoint(e);
        if(noNeighbor == it2)
            return op.InsertPointOnBoundaryEdge(pos, it1, e, outs);
        else return op.InsertPointOnSharedEdge(pos, it1, it2, outs);
    }

    bool isNeedleTriangle(TriIdx it) const
    {
        float_t threshold = 1e3;
        auto angles = Base::Utility::GetInteriorAngles(tri, it);
        std::sort(angles.begin(), angles.end(), std::less<float_t>());
        return (angles[2] / angles[0] > threshold) && (angles[1] / angles[0] > threshold);
    }

    void RemoveNeedleTriangle(TriIdx it)
    {
        auto e = Base::Utility::GetMinLenEdge(tri, it);
        op.TryDecayOneEdge(e);
    }

private:
    VerIdxSet m_circumcenters;
    mutable EncroachedEdgeQueue m_queueE;
    mutable SkinnyTriangleQueue m_queueT;
};

template <typename num_type>
class RuppertRefinement2D : public DelaunayRefinement2D<num_type>
{
    using VerIdxMap = std::unordered_map<VerIdx, VerIdx>;
    using TriIdxMap = std::unordered_map<TriIdx, TriIdx>;
    using Base = DelaunayRefinement2D<num_type>;
    using float_t = typename Base::float_t;
    using EdgeSet = typename Base::EdgeSet;
    using PointF = typename Base::PointF;
    using Point = typename Base::Point;
    using Base::state;
    using Base::paras;
    using Base::tri;
    using Base::op;
public:
    explicit RuppertRefinement2D(Triangulation<Point> & t)
     : DelaunayRefinement2D<num_type>(t)
    {
        Base::SetParas(LimitAlpha(), 0, 0);
    }

    float_t LimitAlpha() const { return math::Rad(20); };

    void Refine(index_t maxStep)
    {
        while(maxStep && (state.curr.fe || state.curr.ft)){
            if(state.curr.fe){
                SplitEncroachedEdge(tri.VerIdxPoint(state.curr.ive), state.curr.e);
            }
            else{
                IndexEdge e(noVertex, noVertex);
                PointF cc = Base::Utility::GetCircumCenter(tri, state.curr.its);
                if(FindFirstEncroachedEdge(cc.template Cast<num_type>(), e)){
                    SplitEncroachedEdge(cc.template Cast<num_type>(), e);
                }
                else{
                    const auto & triangle = tri.triangles[state.curr.its];
                    auto loc = Base::Utility::GetPointTriangleLocation(tri, state.curr.its, cc);
                    if(PointTriangleLocation::Inside == loc) InsertCircumCenter(cc.template Cast<num_type>());
                    else if(PointTriangleLocation::Outside == loc) InsertCircumCenter(cc.template Cast<num_type>());
                    else SplitEncroachedEdge(cc.template Cast<num_type>(), triangle.Edge(static_cast<int>(loc) - 1));
                }
            }
            UpdateState();
            maxStep--;
        }
    }
private:
    void UpdateState()
    {
        state.curr.fe = FindFirstEncroachedEdge(state.curr.e, state.curr.ive);
        state.curr.ft = Base::FindFirstSkinnyTriangle(state.curr.its);
    }

    void SplitEncroachedEdge(const Point & p, const IndexEdge & e)
    {
        auto [e1, e2] = SplitEdge(e);
        if(Base::stillExist(e1) && isPointEncroached(p, e1))
            SplitEncroachedEdge(p, e1);
        if(Base::stillExist(e2) && isPointEncroached(p, e2))
            SplitEncroachedEdge(p, e2);
    }

    Point GetEdgeSplitPosition(const IndexEdge & e)
    {
        std::list<IndexEdge> edges;
        const auto & p1 = tri.VerIdxPoint(e.v1());
        const auto & p2 = tri.VerIdxPoint(e.v2());
        edges = tri.ConcentricFixedEdges(e.v1());
        if(edges.size() > 1){
            for(const auto & edge : edges){
                if(edge == e) continue;
                const auto & p = edge.v1() == e.v1() ? tri.VerIdxPoint(edge.v2()) : tri.VerIdxPoint(edge.v1());
                if(Angle(p, p1, p2) < math::pi_half)
                    return Base::SplitEdgeAtCircularShell(e);
            }
        }

        edges = tri.ConcentricFixedEdges(e.v2());
        if(edges.size() > 1){
            for(const auto & edge : edges){
                if(edge == e) continue;
                const auto & p = edge.v1() == e.v2() ? tri.VerIdxPoint(edge.v2()) : tri.VerIdxPoint(edge.v1());
                if(Angle(p, p2, p1) < math::pi_half)
                    return Base::SplitEdgeAtCircularShell(e);
            }
        }

        return Base::SplitEdgeAtMidPoint(e);
    }

    std::pair<IndexEdge, IndexEdge> SplitEdge(const IndexEdge & e)
    {
        auto split = Base::SplitEdgeAtMidPoint(e);//test, bwu
        //auto split = GetEdgeSplitPosition(e);

        VerIdx iv;
        std::stack<TriIdx> tris;
        auto [it, itOp] = tri.GetTriangles(e);
        if(itOp == noNeighbor)
            iv = op.InsertPointOnBoundaryEdge(split, it, e, tris);
        else iv = op.InsertPointOnSharedEdge(split, it, itOp, tris);

        auto sub = std::make_pair(IndexEdge(e.v1(), iv), IndexEdge(iv, e.v2()));
        if(tri.fixedEdges.count(e)){
            tri.fixedEdges.erase(e);
            tri.fixedEdges.insert(sub.first);
            tri.fixedEdges.insert(sub.second);
        }

        op.Delaunay(tris, iv);
        return sub;
    }

    void InsertCircumCenter(const Point & c)
    {
        VerIdx iv;
        std::stack<TriIdx> tris;
        auto [it1, it2] = op.TraversalTriangleAt(c);
        if(noNeighbor == it2){
            auto loc = Base::Utility::GetPointTriangleLocation(tri, it1, c);
            GENERIC_ASSERT(loc != PointTriangleLocation::Outside)
            if(PointTriangleLocation::Inside == loc)
                iv = op.InsertPointInTriangle(c, it1, tris);
            else {
                auto edge = tri.triangles[it1].Edge(static_cast<int>(loc) - 1);
                iv = op.InsertPointOnBoundaryEdge(c, it1, edge, tris);
            }
        }
        else iv = op.InsertPointOnSharedEdge(c, it1, it2, tris);
        op.Delaunay(tris, iv);
    }

    bool hasSeditiousEdge(const TriIdx & it) const
    {
        for(auto i = 0; i < 3; ++i){
            auto j = (i + 1) % 3;
            auto k = (j + 1) % 3;
            const auto & p1 = tri.TriIdxPoint(it, i);
            const auto & p2 = tri.TriIdxPoint(it, j);
            const auto & p3 = tri.TriIdxPoint(it, k);
            if(Angle(p1, p2, p3) < math::pi / 3){
                size_t l = 0;
                float_t d1 = Distance(p1, p2);
                float_t d2 = Distance(p2, p3);
                float_t d3 = Distance(p3, p1);
                for(; l < sizeof(num_type) * 8; ++l){
                    if(std::pow(2, l) < d1 && d1 < std::pow(2, l + 1)) break;
                }
                if(std::pow(2, l) < d2 && d2 < std::pow(2, l + 1)){
                    //std::cout << "has seditious edge!" << std::endl;//test, bwu
                    return true;
                }
            }
        }
        return false;
    }

    bool FindFirstEncroachedEdge(IndexEdge & edge, VerIdx & ive) const
    {
        //todo, refactory
        for(const auto & e : tri.fixedEdges){
            const auto & tris = tri.GetTriangles(e);
            if(tris.first == noNeighbor){
                continue;//test, bwu
            }
            const auto & t1 = tri.triangles.at(tris.first);
            const auto & v1 = t1.EgOpVe(e);
            const auto & p1 = tri.VerIdxPoint(v1);
            if(isPointEncroached(p1, e)){
                edge = e; ive = v1;
                return true;
            }
            if(tris.second != noNeighbor){
                const auto & t2 = tri.triangles.at(tris.second);
                const auto & v2 = t2.EgOpVe(e);
                const auto & p2 = tri.VerIdxPoint(v2);
                if(isPointEncroached(p2, e)){
                    edge = e; ive = v2;
                    return true;
                }
            }
        }
        return false;
    }

    bool FindFirstEncroachedEdge(const Point & pos, IndexEdge & edge) const
    {
        for(const auto & triangle : tri.triangles){
            for(index_t ie = 0; ie < 3; ++ie){
                edge = triangle.Edge(ie);
                if(isPointEncroached(pos, edge))
                    return true;
            }
        }
        return false;
    }

    bool isPointEncroached(const Point & p, const IndexEdge & e) const
    {
        const auto & ps = tri.VerIdxPoint(e.v1());
        const auto & pe = tri.VerIdxPoint(e.v2());

        auto mid = (ps + pe) * 0.5;
        auto diaSq = (ps - pe).NormSquare();
        if(diaSq < paras.minEdgeLenSq) return false;////test, bwu

        auto disSq = (p - mid).NormSquare();
        if(disSq * 4 < diaSq) return true;
        return false;
    }
};

template <typename num_type>
class ChewSecondRefinement2D : public DelaunayRefinement2D<num_type>
{
    using VerIdxMap = std::unordered_map<VerIdx, VerIdx>;
    using TriIdxMap = std::unordered_map<TriIdx, TriIdx>;
    using Base = DelaunayRefinement2D<num_type>;
    using float_t = typename Base::float_t;
    using EdgeSet = typename Base::EdgeSet;
    using PointF = typename Base::PointF;
    using Point = typename Base::Point;
    using Base::state;
    using Base::paras;
    using Base::tri;
    using Base::op;
    
    struct SkinnyTriangle
    {
        TriIdx index;
        float_t weight;
    };

    struct SkinnyTriangleCmp
    {
        bool operator() (const SkinnyTriangle & t1, const SkinnyTriangle & t2) const noexcept
        {
            return t1.weight < t2.weight;
        }
    };
    using SkinnyTriangleQueue = std::priority_queue<SkinnyTriangle, std::vector<SkinnyTriangle>, SkinnyTriangleCmp>;
public:
    explicit ChewSecondRefinement2D(Triangulation<Point> & t)
     : DelaunayRefinement2D<num_type>(t)
    {
        Base::SetParas(LimitAlpha(), 0, 0);
    }

    float_t LimitAlpha() const { return math::Rad(27.5); }

    void Refine(index_t maxStep)
    {
        index_t step = 0;
        while(maxStep && state.curr.ft){
            if(state.curr.fe){
                SliptEncroachedEdgeAndRemoveCircumCenters(state.curr.e);
            }
            else {
                InsertCircumCenter(state.curr.its);
            }
            UpdateState();
            maxStep--;
        }
    }
private:
    void UpdateState()
    {
        state.curr.ft = UpdateFirstSkinnyTriangle(state.curr.its);
        state.curr.fe = FindEncroachedFixedEdge(state.curr.its, state.curr.e, state.curr.onEdge);
    }

    bool UpdateFirstSkinnyTriangle(TriIdx & it)
    {
        EnqueueSkinnyTriangle();
        while(!m_queueT.empty()){
            it = m_queueT.top().index;
            m_queueT.pop();
            if(!op.isTriangleRemoved(it)){
                if(Base::isSkinnyTriangle(it))
                    return true;
            }
            EnqueueSkinnyTriangle();
        }
        return false;
    }

    void EnqueueSkinnyTriangle()
    {
        if(!m_queueT.empty()) return;
        for(auto it = 0; it < tri.triangles.size(); ++it){
            if(op.isTriangleRemoved(it)) continue;
            if(Base::isSkinnyTriangle(it))
                m_queueT.push({it, WT(it)});
        }
    }

    void SliptEncroachedEdgeAndRemoveCircumCenters(const IndexEdge & e)
    {
        SplitEncroachedEdge(e);
        if(!state.curr.onEdge)
            RemoveEncroachedCircumCenters(e);
    }

    VerIdx SplitEncroachedEdge(const IndexEdge & e)
    {
        auto split = Base::SplitEdgeAtMidPoint(e);

        VerIdx iv;
        std::stack<TriIdx> tris;
        auto [it, itOp] = tri.GetTriangles(e);
        if(itOp == noNeighbor)
            iv = op.InsertPointOnBoundaryEdge(split, it, e, tris);
        else iv = op.InsertPointOnSharedEdge(split, it, itOp, tris);
        op.Delaunay(tris, iv);

        if(tri.fixedEdges.count(e)){
            tri.fixedEdges.erase(e);
            tri.fixedEdges.insert(IndexEdge(e.v1(), iv));
            tri.fixedEdges.insert(IndexEdge(iv, e.v2()));
        }
        
        return iv;
    }

    void RemoveEncroachedCircumCenters(const IndexEdge & e)
    {
        EdgeSet free;//test, bwu
        // std::vector<VerIdx> free;//some problem here
        std::list<VerIdx> encroached;
        Base::FindEncroachedCircumCenters(m_circumcenters, e, encroached);
        for(auto iv : encroached){
            Base::RemoveOneVertexFromTriangulation(iv, free);
            Base::ReTriangulate(free);
            m_circumcenters.erase(iv);
        }
    }

    void InsertCircumCenter(const TriIdx & it)
    {
        PointF cc;
        if(Base::isBadlySkinnyTriangle(it))
            cc = Base::GetOffCenterPoint(it);
        else cc = Base::GetCircumCenterPoint(it);

        auto loc = Base::Utility::GetPointTriangleLocation(tri, it, cc);
        if(loc == PointTriangleLocation::Inside && Base::isTinySizedTriangle(it)) return;

        if(PointTriangleLocation::Inside != loc &&
            PointTriangleLocation::Outside != loc){
            const auto & triangle = tri.triangles[it];
            SplitEncroachedEdge(triangle.Edge(static_cast<int>(loc) - 1));
        }
        else m_circumcenters.insert(Base::InsertOnePointInTriangulation(cc.template Cast<typename Point::coor_t>()));
    }

    bool FindEncroachedEdge(const TriIdx & it, IndexEdge & e, bool & onEdge) const
    {
        onEdge = false;
        if(state.curr.ft){
            PointF cc = Base::GetCircumCenterPoint(it);
            auto loc = Base::Utility::GetPointTriangleLocation(tri, it, cc);
            const auto & triangle = tri.triangles[state.curr.its];
            if(PointTriangleLocation::Outside == loc){
                auto tc = Base::Utility::GetCenter(tri, it);
                auto s1 = segment_type<typename Base::PointF>(cc, tc);
                for(const auto & edge : tri.fixedEdges){
                    auto s2 = Base::Utility::GetSegment(tri, edge);
                    if(boost::polygon::intersects(s1, s2)){
                        e = edge;
                        return true;
                    }
                }
            }
            else if(PointTriangleLocation::Inside != loc){
                e = triangle.Edge(static_cast<int>(loc) - 1);
                onEdge = true;
                return true;
            }
        }
        return false;  
    }

    bool isEdgeVisiable2Triangle(const TriIdx & it, const IndexEdge & e) const
    {
        const auto & triangle = tri.triangles[it];
        for(VerIdx iv : triangle.vertices){
            if(iv == e.v1() || iv == e.v2())
                return true;
        }
        return false;
    }

    bool FindEncroachedFixedEdge(const TriIdx & it, IndexEdge & e, bool & onEdge) const
    {
        onEdge = false;
        if(state.curr.ft){
            PointF cc;
            if(Base::isBadlySkinnyTriangle(it))
                cc = Base::GetOffCenterPoint(it);
            else cc = Base::GetCircumCenterPoint(it);
            auto loc = Base::Utility::GetPointTriangleLocation(tri, it, cc);
            const auto & triangle = tri.triangles[state.curr.its];
            if(PointTriangleLocation::Outside == loc){
                auto tc = Base::Utility::GetCenter(tri, it);
                auto s1 = segment_type<PointF>(cc, tc);
                for(const auto & edge : tri.fixedEdges){
                    auto s2 = Base::Utility::GetSegment(tri, edge);
                    if(Intersects(s1, s2)){
                        if(Contains(s2, cc)){
                            onEdge = true;
                        }
                        e = edge;
                        return true;
                    }
                }
            }
            else if(PointTriangleLocation::Inside != loc){
                e = triangle.Edge(static_cast<int>(loc) - 1);
                if(tri.fixedEdges.count(e)){
                    onEdge = true;
                    return true;
                }
            }
        }
        return false;
    }

private:    
    void ReallocateTriangulation()
    {
        VerIdxMap verIdxMap;    
        TriIdxMap triIdxMap;
        if(op.ReallocateTriangulation(verIdxMap, triIdxMap)){
            //remap circumcenters
            VerIdxSet circumcenters;
            for(auto & iv : m_circumcenters)
                circumcenters.insert(verIdxMap.at(iv));
            std::swap(m_circumcenters, circumcenters);

            // remap skinny triangle and encroached edge
           if(state.curr.ft)
               state.curr.its = triIdxMap.at(state.curr.its);
           if(state.curr.fe){
               VerIdx v1 = verIdxMap.at(state.curr.e.v1());
               VerIdx v2 = verIdxMap.at(state.curr.e.v2());
               state.curr.e = IndexEdge(v1, v2);
           }
        }
    }

    float_t WT(TriIdx it) const
    {
        // return Base::Utility::CircumRadius2ShortestEdgeRatioSq(tri, it);//test, bwu
        return Base::Utility::GetTriangleArea(tri, it);
    }

private:
    VerIdxSet m_circumcenters;
    SkinnyTriangleQueue m_queueT;
};

}//namespace tri
}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_TRI_TRIANGULATIONREFINEMENT_HPP
