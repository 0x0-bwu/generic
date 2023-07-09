/**
 * @file Clipper.hpp
 * @author bwu
 * @brief Modified version of clipper library, origin: http://www.angusj.com/delphi/clipper.php
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include <boost/multiprecision/cpp_int.hpp>
#include "generic/common/Exception.hpp"
#include "Polygon.hpp"
#include <queue>

namespace generic  {
namespace geometry {
namespace clipper  {

using namespace common;
using Int128 = boost::multiprecision::int128_t;

namespace constant {
inline static constexpr int64_t loRange = 0x3FFFFFFF;
inline static constexpr int64_t hiRange = 0x3FFFFFFFFFFFFFFFLL;
inline static constexpr int unassigned = -1;  //edge not currently 'owning' a solution
inline static constexpr int skip = -2;        //edge that would otherwise close a path
inline static constexpr double horizontal = -1.0e40;
}//namespace constant

inline Int128 Int128Mul (int64_t lhs, int64_t rhs)
{
    Int128 res = lhs; res *= rhs;
    return res;
};

enum class EndType  { ClosedPolygon, ClosedLine, OpenButt, OpenSquare, OpenRound };
enum class ClipType { Intersection, Union, Difference, Xor };
enum class EdgeSide { Left = 1, Right = 2 };
enum class JoinType { Square, Round, Miter };
enum class PolyType { Subject, Clip };
enum class Direction { RightToLeft, LeftToRight };
enum class InitOptions { ReverseSolution = 1, StrictlySimple = 2, PreserveCollinear = 4};
enum class PolyFillType { EvenOdd, NonZero, Positive, Negative };

template <typename num_type>
using Point = Point2D<num_type>;

template <typename num_type>
using Path = Polyline2D<num_type>;

template <typename num_type>
using Paths = std::vector<Path<num_type> >;

template <typename num_type>
class PolyNode;

template <typename num_type>
using PolyNodes = std::vector<PolyNode<num_type> * >;

template <typename num_type>
class ClipperBase;

template <typename num_type>
class Clipper;

template <typename num_type>
class PolyNode
{
    friend class Clipper<num_type>;
public:
    PolyNode() = default;
    virtual ~PolyNode() {};
    Path<num_type> contour;
    PolyNode<num_type> * parent = nullptr;
    PolyNodes<num_type> children;
    PolyNode<num_type> * GetNext() const;
    bool isHole() const;
    bool isOpen() const;
    size_t ChildCount() const;
private:
    //PolyNode& operator =(PolyNode& other);
    size_t m_index = 0; //node index in Parent.Childs
    bool m_isOpen = false;
    JoinType m_joinType;
    EndType m_endType;
    PolyNode<num_type> * GetNextSiblingUp() const;
    void AddChild(PolyNode & child);
};

template <typename num_type>
class PolyTree : public PolyNode<num_type>
{
    friend class Clipper<num_type>;
public:
    ~PolyTree();
    PolyNode<num_type> * GetFirst() const;
    void Clear();
    size_t Total() const;
private:
    PolyNodes<num_type> m_allNodes;
};

template <typename num_type>
inline PolyNode<num_type> * PolyNode<num_type>::GetNext() const
{
    if(!children.empty()) return children.front();
    else return GetNextSiblingUp();
}

template <typename num_type>
inline PolyNode<num_type> * PolyNode<num_type>::GetNextSiblingUp() const
{
    if(!parent) return nullptr;
    else if(m_index == parent->ChildCount() - 1)
        return parent->GetNextSiblingUp();
    else return parent->children[m_index + 1];
}

template <typename num_type>
inline bool PolyNode<num_type>::isHole() const
{
  bool result = true;
  auto * node = parent;
  while(node) {
      result = !result;
      node = node->parent;
  }
  return result;
}

template <typename num_type>
inline bool PolyNode<num_type>::isOpen() const
{
    return m_isOpen;
}

template <typename num_type>
inline size_t PolyNode<num_type>::ChildCount() const
{
    return children.size();
}

template <typename num_type>
inline void PolyNode<num_type>::AddChild(PolyNode<num_type> & child)
{
    auto index = ChildCount();
    children.push_back(&child);
    child.parent = this;
    child.m_index = index;
}

template <typename num_type>
inline PolyTree<num_type>::~PolyTree()
{
    Clear();
}

template <typename num_type>
inline void PolyTree<num_type>::Clear()
{
    for(auto & node : m_allNodes)
        delete node;

    m_allNodes.clear();
    PolyNode<num_type>::children.clear();
}

template <typename num_type>
inline PolyNode<num_type> * PolyTree<num_type>::GetFirst() const
{
    if(!PolyNode<num_type>::children.empty())
        return PolyNode<num_type>::children.front();
    return nullptr;
}

template <typename num_type>
inline size_t PolyTree<num_type>::Total() const
{
    auto result = m_allNodes.size();
    if(result > 0 &&
        PolyNode<num_type>::children.front() != m_allNodes.front())
        result--;
    return result;
}

template <typename num_type>
struct OutPt
{
    int idx;
    Point<num_type> pt;
    OutPt * next = nullptr;
    OutPt * prev = nullptr;
};

//OutRec: contains a path in the clipping solution. Edges in the AEL will
//carry a pointer to an OutRec when they are part of the clipping solution.
template <typename num_type>
struct OutRec {
  int idx;
  bool isHole = false;
  bool isOpen = false;
  OutRec<num_type> * firstLeft = nullptr;  //see comments in clipper.pas
  PolyNode<num_type> * polyNd = nullptr;
  OutPt<num_type> * pts = nullptr;
  OutPt<num_type> * bottomPt = nullptr;
};

template <typename num_type>
struct Join {
    OutPt<num_type> * outPt1 = nullptr;
    OutPt<num_type> * outPt2 = nullptr;
    Point<num_type> offPt;
};

template <typename num_type>
struct TEdge
{
    using float_t = float_type<num_type>;
    Point<num_type> bot;
    Point<num_type> curr; //current (updated for every new scanbeam)
    Point<num_type> top;
    float_t dx;
    PolyType polyTyp;
    EdgeSide side; //side only refers to current side of solution poly
    int windDelta; //1 or -1 depending on winding direction
    int windCnt;
    int windCnt2; //winding count of the opposite polytype
    int outIdx;
    TEdge<num_type> * next = nullptr;
    TEdge<num_type> * prev = nullptr;
    TEdge<num_type> * nextInLML = nullptr;
    TEdge<num_type> * nextInAEL = nullptr;
    TEdge<num_type> * prevInAEL = nullptr;
    TEdge<num_type> * nextInSEL = nullptr;
    TEdge<num_type> * prevInSEL = nullptr;
};

template <typename num_type>
struct IntersectNode
{
    TEdge<num_type> * edge1 = nullptr;
    TEdge<num_type> * edge2 = nullptr;
    Point<num_type> point;
};

template <typename num_type>
using IntersectList = std::vector<IntersectNode<num_type> * >;

template <typename num_type>
struct LocalMinimum {
    num_type y;
    TEdge<num_type> * boundL = nullptr;
    TEdge<num_type> * boundR = nullptr;
};

template <typename num_type>
inline void InitEdge(TEdge<num_type> * e, TEdge<num_type> * eNext, TEdge<num_type> * ePrev, const Point<num_type> & p)
{
    std::memset(e, 0, sizeof(TEdge<num_type>));
    e->next = eNext;
    e->prev = ePrev;
    e->curr = p;
    e->outIdx = constant::unassigned;
}

template <typename num_type>
inline void SetDx(TEdge<num_type> & e)
{
    using float_t = float_type<num_type>;
    num_type dy = (e.top[1] - e.bot[1]);
    if(dy == 0) e.dx = constant::horizontal;
    else e.dx = static_cast<float_t>(e.top[0] - e.bot[0]) / dy;
}

template <typename num_type>
inline void InitEdge2(TEdge<num_type> & e, PolyType polyTyp)
{
    if(math::GE(e.curr[1], e.next->curr[1])) {
        e.bot = e.curr;
        e.top = e.next->curr;
    }
    else {
        e.top = e.curr;
        e.bot = e.next->curr;
    }
    SetDx(e);
    e.polyTyp = polyTyp;
}

template <typename num_type>
inline void ReverseHorizontal(TEdge<num_type> & e)
{
    //swap horizontal edges' Top and Bottom x's so they follow the natural
    //progression of the bounds - ie so their xbots will align with the
    //adjoining lower edge. [Helpful in the ProcessHorizontal() method.]
    std::swap(e.top[0], e.bot[0]);
}

template <typename num_type>
inline bool isHorizontal(TEdge<num_type> & e)
{
    using float_t = float_type<num_type>;
    return math::EQ<float_t>(e.dx, constant::horizontal);
}

template <typename num_type>
inline bool HorzSegmentsOverlap(num_type s1a, num_type s1b, num_type s2a, num_type s2b)
{
    if(math::GT(s1a, s1b)) std::swap(s1a, s1b);
    if(math::GT(s2a, s2b)) std::swap(s2a, s2b);
    return math::LT(s1a, s2b) && math::LT(s2a, s1b);
}

template <typename num_type>
inline TEdge<num_type> * FindNextLocMin(TEdge<num_type> * e)
{
    for(;;) {
        while (e->bot != e->prev->bot || e->curr == e->top) e = e->next;
        if (!isHorizontal(*e) && !isHorizontal(*e->prev)) break;
        while(isHorizontal(*e->prev)) e = e->prev;
        auto * e2 = e;
        while(isHorizontal(*e)) e = e->next;
        if(math::EQ(e->top[1], e->prev->bot[1])) continue; //ie just an intermediate horz.
        if(math::LT(e2->prev->bot[0], e->bot[0])) e = e2;
        break;
    }
    return e;
}

template <typename num_type>
inline TEdge<num_type> * RemoveEdge(TEdge<num_type> * e)
{
    //removes e from double_linked_list (but without removing from memory)
    e->prev->next = e->next;
    e->next->prev = e->prev;
    auto * result = e->next;
    e->prev = nullptr; //flag as removed (see ClipperBase.Clear)
    return result;
}

template <typename num_type>
inline num_type TopX(TEdge<num_type> & edge, const num_type currentY)
{
  return math::EQ(currentY, edge.top[1]) ?
    edge.top[0] : edge.bot[0] + std::round(edge.dx * (currentY - edge.bot[1]));
}

template <typename num_type>
inline void IntersectPoint(TEdge<num_type> & edge1, TEdge<num_type> & edge2, Point<num_type> & ip)
{
    using float_t = float_type<num_type>;
    float_t b1, b2;
    if(math::EQ(edge1.dx, edge2.dx)) {
        ip[1] = edge1.curr[1];
        ip[0] = TopX(edge1, ip[1]);
        return;
    }
    else if(math::EQ<float_t>(edge1.dx, 0)) {
        ip[0] = edge1.bot[0];
        if(isHorizontal(edge2))
            ip[1] = edge2.bot[1];
        else {
            b2 = edge2.bot[1] - (edge2.bot[0] / edge2.dx);
            ip[1] = std::round(ip[0] / edge2.dx + b2);
        }
    }
    else if(math::EQ<float_t>(edge2.dx, 0)) {
        ip[0] = edge2.bot[0];
        if(isHorizontal(edge1))
            ip[1] = edge1.bot[1];
        else {
            b1 = edge1.bot[1] - (edge1.bot[0] / edge1.dx);
            ip[1] = std::round(ip[0] / edge1.dx + b1);
        }
    }
    else {
        b1 = edge1.bot[0] - edge1.bot[1] * edge1.dx;
        b2 = edge2.bot[0] - edge2.bot[1] * edge2.dx;
        float_t q = (b2 - b1) / (edge1.dx - edge2.dx);
        ip[1] = std::round(q);
        if(std::fabs(edge1.dx) < std::fabs(edge2.dx))
            ip[0] = std::round(edge1.dx * q + b1);
        else ip[0] = std::round(edge2.dx * q + b2);
    }

    if(math::LT(ip[1], edge1.top[1]) || math::LT(ip[1], edge2.top[1])) {
        if(math::GT(edge1.top[1], edge2.top[1]))
            ip[1] = edge1.top[1];
        else ip[1] = edge2.top[1];
        if(math::LT(std::fabs(edge1.dx), std::fabs(edge2.dx)))
            ip[0] = TopX(edge1, ip[1]);
        else ip[0] = TopX(edge2, ip[1]);
    }
    //finally, don't allow 'ip' to be BELOW curr[1] (ie bottom of scanbeam) ...
    if(math::GT(ip[1], edge1.curr[1])) {
        ip[1] = edge1.curr[1];
        //use the more vertical edge to derive X ...
        if(math::GT(std::fabs(edge1.dx), std::fabs(edge2.dx)))
            ip[0] = TopX(edge2, ip[1]);
        else ip[0] = TopX(edge1, ip[1]);
    }
}

template <typename num_type>
inline bool E2InsertsBeforeE1(TEdge<num_type> & e1, TEdge<num_type> & e2)
{
  if(e2.curr[0] == e1.curr[0]) {
    if(e2.top[1] > e1.top[1])
        return math::LT(e2.top[0], TopX(e1, e2.top[1]));
    else return math::GT(e1.top[0], TopX(e2, e1.top[1]));
  }
  else return math::LT(e2.curr[0], e1.curr[0]);
}

template <typename num_type>
inline void RangeTest(const Point<num_type> & p, bool & useFullRange)
{
  if(useFullRange) {
    if (p[0] > constant::hiRange ||  p[1] > constant::hiRange||
       -p[0] > constant::hiRange || -p[1] > constant::hiRange)
      throw Exception("Coordinate outside allowed range");
  }
  else if(p[0] > constant::loRange ||  p[1] > constant::loRange || 
         -p[0] > constant::loRange || -p[1] > constant::loRange) {
    useFullRange = true;
    RangeTest(p, useFullRange);
  }
}

template <typename num_type>
inline bool SlopesEqual(const TEdge<num_type> & e1, const TEdge<num_type> & e2, bool useFullInt64Range)
{
    if constexpr (std::is_same<num_type, int64_t>::value){
        if(useFullInt64Range)
            return Int128Mul(e1.top[1] - e1.bot[1], e2.top[0] - e2.bot[0]) == Int128Mul(e1.top[0] - e1.bot[0], e2.top[1] - e2.bot[1]);
        else
            return math::EQ((e1.top[1] - e1.bot[1]) * (e2.top[0] - e2.bot[0]), (e1.top[0] - e1.bot[0]) * (e2.top[1] - e2.bot[1]));
    }
    else return math::EQ((e1.top[1] - e1.bot[1]) * (e2.top[0] - e2.bot[0]), (e1.top[0] - e1.bot[0]) * (e2.top[1] - e2.bot[1]));
}

template <typename num_type>
inline bool SlopesEqual(const Point<num_type> p1, const Point<num_type> p2, const Point<num_type> p3, bool useFullInt64Range)
{
    if constexpr (std::is_same<num_type, int64_t>::value) {
        if(useFullInt64Range)
            return Int128Mul(p1[1] - p2[1], p2[0] - p3[0]) == Int128Mul(p1[0] - p2[0], p2[1] - p3[1]);
        else return math::EQ((p1[1] - p2[1]) * (p2[0] - p3[0]), (p1[0] - p2[0]) * (p2[1] - p3[1]));
    }
    else return math::EQ((p1[1] - p2[1]) * (p2[0] - p3[0]), (p1[0] - p2[0]) * (p2[1] - p3[1]));
}

template <typename num_type>
inline bool SlopesEqual(const Point<num_type> p1, const Point<num_type> p2, const Point<num_type> p3, const Point<num_type> p4, bool useFullInt64Range)
{
    if constexpr (std::is_same<num_type, int64_t>::value) {
        if(useFullInt64Range)
            return Int128Mul(p1[1] - p2[1], p3[0] - p4[0]) == Int128Mul(p1[0] - p2[0], p3[1] - p4[1]);
        else return math::EQ((p1[1] - p2[1]) * (p3[0] - p4[0]), (p1[0] - p2[0]) * (p3[1] - p4[1]));
    }
    else return math::EQ((p1[1] - p2[1]) * (p3[0] - p4[0]), (p1[0] - p2[0]) * (p3[1] - p4[1]));
}

template <typename num_type>
inline bool Pt2IsBetweenPt1AndPt3(const Point<num_type> & p1, const Point<num_type> & p2, const Point<num_type> & p3)
{
  if((p1 == p3) || (p1 == p2) || (p3 == p2)) return false;
  else if (p1[0] != p3[0]) return (p2[0] > p1[0]) == (p2[0] < p3[0]);
  else return (p2[1] > p1[1]) == (p2[1] < p3[1]);
}

template <typename num_type>
inline void GetHorzDirection(TEdge<num_type> & horzEdge, Direction & dir, num_type & left, num_type & right)
{
    if(math::LT(horzEdge.bot[0], horzEdge.top[0])) {
        left = horzEdge.bot[0];
        right = horzEdge.top[0];
        dir = Direction::LeftToRight;
    }
    else {
        left = horzEdge.top[0];
        right = horzEdge.bot[0];
        dir = Direction::RightToLeft;
    }
}

template <typename num_type>
inline bool isIntermediate(TEdge<num_type> * e, const num_type y)
{
    return math::EQ(e->top[1], y) && e->nextInLML;
}

template <typename num_type>
inline TEdge<num_type> * GetMaximaPair(TEdge<num_type> * e)
{
    if((e->next->top == e->top) && !e->next->nextInLML)
        return e->next;
    else if((e->prev->top == e->top) && !e->prev->nextInLML)
        return e->prev;
    else return nullptr;
}

template <typename num_type>
inline TEdge<num_type> * GetMaximaPairEx(TEdge<num_type> * e)
{
    //as GetMaximaPair() but returns 0 if MaxPair isn't in AEL (unless it's horizontal)
    auto * result = GetMaximaPair(e);
    if(result && (result->outIdx == constant::skip ||
        (result->nextInAEL == result->prevInAEL && !isHorizontal(*result)))) return nullptr;
    return result;
}

template <typename num_type>
inline TEdge<num_type> * GetNextInAEL(TEdge<num_type> * e, Direction dir)
{
    return dir == Direction::LeftToRight ? e->nextInAEL : e->prevInAEL;
}

template <typename num_type>
inline float_type<num_type> GetDx(const Point<num_type> p1, const Point<num_type> p2)
{
  return math::EQ(p1[1], p2[1]) ? constant::horizontal : float_type<num_type>(p2[0] - p1[0]) / (p2[1] - p1[1]);
}

template <typename num_type>
inline bool FirstIsBottomPt(const OutPt<num_type> * btmPt1, const OutPt<num_type> * btmPt2)
{
    using float_t = float_type<num_type>;

    auto * p = btmPt1->prev;
    while(p->pt == btmPt1->pt && p != btmPt1) p = p->prev;
    float_t dx1p = std::fabs(GetDx(btmPt1->pt, p->pt));
    p = btmPt1->next;
    while(p->pt == btmPt1->pt && p != btmPt1) p = p->next;
    float_t dx1n = std::fabs(GetDx(btmPt1->pt, p->pt));

    p = btmPt2->prev;
    while(p->pt == btmPt2->pt && p != btmPt2) p = p->prev;
    float_t dx2p = std::fabs(GetDx(btmPt2->pt, p->pt));
    p = btmPt2->next;
    while(p->pt == btmPt2->pt && p != btmPt2) p = p->next;
    float_t dx2n = std::fabs(GetDx(btmPt2->pt, p->pt));

    if(math::EQ(std::max(dx1p, dx1n), std::max(dx2p, dx2n)) &&
        math::EQ(std::min(dx1p, dx1n), std::min(dx2p, dx2n)))
        return math::isPositive(Area(btmPt1)); //if otherwise identical use orientation
    else return (math::GE(dx1p, dx2p) && math::GE(dx1p, dx2n)) || (math::GE(dx1n, dx2p) && math::GE(dx1n, dx2n));
}

template <typename num_type>
inline OutPt<num_type> * GetBottomPt(OutPt<num_type> * pp)
{
    OutPt<num_type> * dups = nullptr;
    auto * p = pp->next;
    while(p != pp) {
        if(p->pt[1] > pp->pt[1]) {
            pp = p;
            dups = nullptr;
        }
        else if(math::EQ(p->pt[1], pp->pt[1]) && math::LE(p->pt[0], pp->pt[0])) {
            if(math::LT(p->pt[0], pp->pt[0])) {
                dups = nullptr;
                pp = p;
            }
            else {
                if(p->next != pp && p->prev != pp) dups = p;
            }
        }
        p = p->next;
    }

    if(dups) {
        //there appears to be at least 2 vertices at BottomPt so ...
        while(dups != p) {
            if(!FirstIsBottomPt(p, dups)) pp = dups;
            dups = dups->next;
            while(dups->pt != pp->pt) dups = dups->next;
        }
    }
    return pp;
}

template <typename num_type>
inline OutRec<num_type> * GetLowermostRec(OutRec<num_type> * outRec1, OutRec<num_type> * outRec2)
{
    //work out which polygon fragment has the correct hole state ...
    if(!outRec1->bottomPt)
        outRec1->bottomPt = GetBottomPt(outRec1->pts);
    if(!outRec2->bottomPt)
        outRec2->bottomPt = GetBottomPt(outRec2->pts);
    auto * outPt1 = outRec1->bottomPt;
    auto * outPt2 = outRec2->bottomPt;
    if (outPt1->pt[1] > outPt2->pt[1]) return outRec1;
    else if (outPt1->pt[1] < outPt2->pt[1]) return outRec2;
    else if (outPt1->pt[0] < outPt2->pt[0]) return outRec1;
    else if (outPt1->pt[0] > outPt2->pt[0]) return outRec2;
    else if (outPt1->next == outPt1) return outRec2;
    else if (outPt2->next == outPt1) return outRec1;
    else if (FirstIsBottomPt(outPt1, outPt2)) return outRec1;
    else return outRec2;
}

template <typename num_type>
inline bool OutRec1RightOfOutRec2(OutRec<num_type> * outRec1, OutRec<num_type> * outRec2)
{
    do {
        outRec1 = outRec1->firstLeft;
        if(outRec1 == outRec2) return true;
    } while(outRec1);
    return false;
}

template <typename num_type>
inline void ReversePolyPtLinks(OutPt<num_type> * pp)
{
    if(!pp) return;
    OutPt<num_type> *pp1 = nullptr, * pp2 = nullptr;
    pp1 = pp;
    do {
        pp2 = pp1->next;
        pp1->next = pp1->prev;
        pp1->prev = pp2;
        pp1 = pp2;
    } while(pp1 != pp);
}

template <typename num_type>
inline void DisposeOutPts(OutPt<num_type> * & pp)
{
    if(pp == nullptr) return;
    pp->prev->next = nullptr;
    while(pp){
        auto * tmp = pp;
        pp = pp->next;
        delete tmp;
    }
}

template <typename num_type>
inline bool EdgesAdjacent(const IntersectNode<num_type> & node)
{
  return node.edge1->nextInSEL == node.edge2 || node.edge1->prevInSEL == node.edge2;
}

template <typename num_type>
inline bool isMaxima(TEdge<num_type> * e, const num_type y)
{
    return e && math::EQ(e->top[1], y) && !e->nextInLML;
}

template <typename num_type>
inline size_t PointCount(OutPt<num_type> * pts)
{
    if (!pts) return 0;
    size_t ct = 0;
    auto * p = pts;
    do {
        ct++;
        p = p->next;
    }
    while (p != pts);
    return ct;
}

template <typename num_type>
inline float_type<num_type> Area(const OutPt<num_type> * op)
{
    using float_t = float_type<num_type>;
    float_t area = 0;
    const auto * startOp = op;
    if(nullptr == op) return area;
    do {
        area += float_t(op->prev->pt[0] + op->pt[0]) * (op->prev->pt[1] - op->pt[1]);
        op = op->next;
    } while(op != startOp);
    return 0.5 * area;
}

template <typename num_type>
inline float_type<num_type> Area(const OutRec<num_type> & outRec)
{
    return Area(outRec.pts);
}

template <typename num_type>
inline OutPt<num_type> * DupOutPt(OutPt<num_type> * outPt, bool insertAfter)
{
    auto * result = new OutPt<num_type>;
    result->pt = outPt->pt;
    result->idx = outPt->idx;
    if(insertAfter) {
        result->next = outPt->next;
        result->prev = outPt;
        outPt->next->prev = result;
        outPt->next = result;
    }
    else {
        result->prev = outPt->prev;
        result->next = outPt;
        outPt->prev->next = result;
        outPt->prev = result;
    }
    return result;
}

template <typename num_type>
inline void UpdateOutPtIdxs(OutRec<num_type> & outrec)
{
    auto * op = outrec.pts;
    do {
        op->idx = outrec.idx;
        op = op->prev;
    } while(op != outrec.pts);
}

template <typename num_type>
inline int PointInPolygon (const Point<num_type> & pt, OutPt<num_type> * op)
{
    using float_t = float_type<num_type>;
    //returns 0 if false, +1 if true, -1 if pt ON polygon boundary
    int result = 0;
    auto * startOp = op;
    for(;;) {
        if(op->next->pt[1] == pt[1]) {
            if(math::EQ(op->next->pt[0], pt[0]) ||
                (math::EQ(op->pt[1], pt[1]) && math::GT(op->next->pt[0], pt[0]) == math::LT(op->pt[0], pt[0]))) return -1;
        }
        if(math::LT(op->pt[1], pt[1]) != math::LT(op->next->pt[1], pt[1])) {
            if(math::GE(op->pt[0], pt[0])) {
                if(math::GT(op->next->pt[0], pt[0]))
                    result = 1 - result;
                else {
                    float_t d = float_t(op->pt[0] - pt[0]) * (op->next->pt[1] - pt[1]) - float_t(op->next->pt[0] - pt[0]) * (op->pt[1] - pt[1]);
                    if(math::EQ<float_t>(d, 0)) return -1;
                    if(math::isPositive(d) == math::GT(op->next->pt[1], op->pt[1])) result = 1 - result;
                }
            }
            else {
                if(math::GT(op->next->pt[0], pt[0])) {
                    float_t d = float_t(op->pt[0] - pt[0]) * (op->next->pt[1] - pt[1]) - float_t(op->next->pt[0] - pt[0]) * (op->pt[1] - pt[1]);
                    if(math::EQ<float_t>(d, 0)) return -1;
                    if(math::isPositive(d) == math::GT(op->next->pt[1], op->pt[1])) result = 1 - result;
                }
            }
        }
        op = op->next;
        if(startOp == op) break;
    }
    return result;
}

template <typename num_type>
inline bool Poly2ContainsPoly1(OutPt<num_type> * outPt1, OutPt<num_type> * outPt2)
{
    auto * op = outPt1;
    do {
        //nb: PointInPolygon returns 0 if false, +1 if true, -1 if pt on polygon
        int res = PointInPolygon(op->pt, outPt2);
        if (res >= 0) return res > 0;
        op = op->next;
    } while (op != outPt1);
    return true;
}

template <typename num_type>
inline OutRec<num_type> * ParseFirstLeft(OutRec<num_type> * firstLeft)
{
    while(firstLeft && !firstLeft->pts)
        firstLeft = firstLeft->firstLeft;
    return firstLeft;
}

template <typename num_type>
inline bool GetOverlap(const num_type a1, const num_type a2, const num_type b1, const num_type b2, num_type & left, num_type & right)
{
    if(math::LT(a1, a2)) {
        if(math::LT(b1, b2)) { left = std::max(a1, b1); right = std::min(a2, b2); }
        else { left = std::max(a1, b2); right = std::min(a2, b1); }
    }
    else {
        if (math::LT(b1, b2)) { left = std::max(a2, b1); right = std::min(a1, b2); }
        else { left = std::max(a2, b2); right = std::min(a1, b1); }
    }
    return math::LT(left, right);
}

template <typename num_type>
inline bool JoinHorz(OutPt<num_type> * op1, OutPt<num_type> * op1b, OutPt<num_type> * op2, OutPt<num_type> * op2b, const Point<num_type> pt, bool discardLeft)
{
    auto dir1 = math::GT(op1->pt[0], op1b->pt[0]) ? Direction::RightToLeft : Direction::LeftToRight;
    auto dir2 = math::GT(op2->pt[0], op2b->pt[0]) ? Direction::RightToLeft : Direction::LeftToRight;
    if(dir1 == dir2) return false;

    //When DiscardLeft, we want Op1b to be on the Left of Op1, otherwise we
    //want Op1b to be on the Right. (And likewise with Op2 and Op2b.)
    //So, to facilitate this while inserting Op1b and Op2b ...
    //when DiscardLeft, make sure we're AT or RIGHT of Pt before adding Op1b,
    //otherwise make sure we're AT or LEFT of Pt. (Likewise with Op2b.)
    if(dir1 == Direction::LeftToRight) {
        while(math::LE(op1->next->pt[0], pt[0]) &&
                math::GE(op1->next->pt[0], op1->pt[0]) && math::EQ(op1->next->pt[1], pt[1]))
            op1 = op1->next;
        if(discardLeft && math::NE(op1->pt[0], pt[0])) op1 = op1->next;
        op1b = DupOutPt(op1, !discardLeft);
        if(op1b->pt != pt) {
            op1 = op1b;
            op1->pt = pt;
            op1b = DupOutPt(op1, !discardLeft);
        }
    }
    else {
        while (math::GE(op1->next->pt[0], pt[0]) &&
                math::LE(op1->next->pt[0], op1->pt[0]) && math::EQ(op1->next->pt[1], pt[1]))
            op1 = op1->next;
        if(!discardLeft && math::NE(op1->pt[0], pt[0])) op1 = op1->next;
        op1b = DupOutPt(op1, discardLeft);
        if(op1b->pt != pt) {
            op1 = op1b;
            op1->pt = pt;
            op1b = DupOutPt(op1, discardLeft);
        }
    }

    if(dir2 == Direction::LeftToRight) {
        while(math::LE(op2->next->pt[0], pt[0]) &&
                math::GE(op2->next->pt[0], op2->pt[0]) && math::EQ(op2->next->pt[1], pt[1]))
            op2 = op2->next;
        if(discardLeft && math::NE(op2->pt[0], pt[0])) op2 = op2->next;
        op2b = DupOutPt(op2, !discardLeft);
        if(op2b->pt != pt) {
            op2 = op2b;
            op2->pt = pt;
            op2b = DupOutPt(op2, !discardLeft);
        };
    }
    else {
        while(math::GE(op2->next->pt[0], pt[0]) &&
                math::LE(op2->next->pt[0], op2->pt[0]) && math::EQ(op2->next->pt[1], pt[1]))
            op2 = op2->next;
        if(!discardLeft && math::NE(op2->pt[0], pt[0])) op2 = op2->next;
        op2b = DupOutPt(op2, discardLeft);
        if(op2b->pt != pt) {
            op2 = op2b;
            op2->pt = pt;
            op2b = DupOutPt(op2, discardLeft);
        };
    };

    if((dir1 == Direction::LeftToRight) == discardLeft) {
        op1->prev = op2;
        op2->next = op1;
        op1b->next = op2b;
        op2b->prev = op1b;
    }
    else {
        op1->next = op2;
        op2->prev = op1;
        op1b->prev = op2b;
        op2b->next = op1b;
    }
    return true;
}

template <typename num_type>
inline void SwapSides(TEdge<num_type> & edge1, TEdge<num_type> & edge2)
{
    std::swap(edge1.side, edge2.side);
}

template <typename num_type>
inline void SwapPolyIndexes(TEdge<num_type> & edge1, TEdge<num_type> & edge2)
{
    std::swap(edge1.outIdx, edge2.outIdx);
}

template <typename num_type>
using EdgeList = std::vector<TEdge<num_type> * >;

template <typename num_type>
using PolyOutList = std::vector<OutRec<num_type> * >;

template <typename num_type>
using JoinList = std::vector<Join<num_type> * >;

template <typename num_type>
class ClipperBase
{
protected:
    using ScanbeamList = std::priority_queue<num_type> ;
    using MinimaList = std::vector<LocalMinimum<num_type> > ;
    using MinimaListIter = typename MinimaList::iterator;

public:
    ClipperBase();
    virtual ~ClipperBase();
    virtual bool AddPath(const Path<num_type> & path, PolyType polyTyp, bool closed);
    bool AddPaths(const Paths<num_type> & paths, PolyType polyTyp, bool closed);
    virtual void Clear();
    // IntRect GetBounds();
    // bool PreserveCollinear() {return m_PreserveCollinear;};
    // void PreserveCollinear(bool value) {m_PreserveCollinear = value;};
protected:
    void DisposeLocalMinimaList();
    // TEdge* AddBoundsToLML(TEdge *e, bool IsClosed);
    virtual void Reset();
    TEdge<num_type> * ProcessBound(TEdge<num_type> * e, bool isClockwise);
    void InsertScanbeam(num_type y);
    bool PopScanbeam(num_type & y);
    bool LocalMinimaPending();
    bool PopLocalMinima(num_type y, const LocalMinimum<num_type> * & locMin);
    OutRec<num_type> * CreateOutRec();
    void DisposeAllOutRecs();
    void DisposeOutRec(size_t index);
    void SwapPositionsInAEL(TEdge<num_type> * edge1, TEdge<num_type> * edge2);
    void DeleteFromAEL(TEdge<num_type> * e);
    void UpdateEdgeIntoAEL(TEdge<num_type> * & e);

protected:
    MinimaListIter m_currentLM;
    MinimaList m_minimaList;
    bool m_useFullRange;
    EdgeList<num_type> m_edges;
    bool m_preserveCollinear;
    bool m_hasOpenPaths;
    PolyOutList<num_type> m_polyOuts;
    TEdge<num_type> * m_activeEdges;
    ScanbeamList m_scanbeam;
};

template <typename num_type>
inline ClipperBase<num_type>::ClipperBase()
{
    m_currentLM = m_minimaList.begin(); //begin() == end() here
    m_useFullRange = false;
}

template <typename num_type>
inline ClipperBase<num_type>::~ClipperBase()
{
    Clear();
}

template <typename num_type>
inline void ClipperBase<num_type>::Clear()
{
    DisposeLocalMinimaList();
    for(size_t i = 0; i < m_edges.size(); ++i) {
        auto * edges = m_edges[i];
        delete [] edges;
    }
    m_edges.clear();
    m_useFullRange = false;
    m_hasOpenPaths = false;
}

template <typename num_type>
inline void ClipperBase<num_type>::DisposeLocalMinimaList()
{
    m_minimaList.clear();
    m_currentLM = m_minimaList.begin();
}

template <typename num_type>
inline void ClipperBase<num_type>::Reset()
{
    m_currentLM = m_minimaList.begin();
    if(m_currentLM == m_minimaList.end()) return; //ie nothing to process

    using LM = LocalMinimum<num_type>;
    std::sort(m_minimaList.begin(), m_minimaList.end(),
            [](const LM & lm1, const LM & lm2){ return math::LT(lm2.y, lm1.y); });

    m_scanbeam = ScanbeamList{}; //clears/resets priority_queue
    //reset all edges ...
    for (const auto & lm : m_minimaList) {
        InsertScanbeam(lm.y);
        auto * e = lm.boundL;
        if(e) {
            e->curr = e->bot;
            e->side = EdgeSide::Left;
            e->outIdx = constant::unassigned;
        }

        e = lm.boundR;
        if(e) {
            e->curr = e->bot;
            e->side = EdgeSide::Right;
            e->outIdx = constant::unassigned;
        }
    }
    m_activeEdges = 0;
    m_currentLM = m_minimaList.begin();
}

template <typename num_type>
inline bool ClipperBase<num_type>::AddPath(const Path<num_type> & path, PolyType polyTyp, bool closed)
{
    if(!closed && polyTyp == PolyType::Clip)
        throw Exception("AddPath: Open paths must be subject.");

    int highI = static_cast<int>(path.size()) - 1;
    if (closed) while (highI > 0 && (path[highI] == path[0])) --highI;
    while (highI > 0 && (path[highI] == path[highI - 1])) --highI;
    if ((closed && highI < 2) || (!closed && highI < 1)) return false;

    //create a new edge array ...
    auto * edges = new TEdge<num_type>[highI + 1];

    bool isFlat = true;
    //1. Basic (first) edge initialization ...
    try
    {
        edges[1].curr = path[1];
        RangeTest(path[0], m_useFullRange);
        RangeTest(path[highI], m_useFullRange);
        InitEdge(&edges[0], &edges[1], &edges[highI], path[0]);
        InitEdge(&edges[highI], &edges[0], &edges[highI-1], path[highI]);
        for(int i = highI - 1; i >= 1; --i) {
            RangeTest(path[i], m_useFullRange);
            InitEdge(&edges[i], &edges[i+1], &edges[i-1], path[i]);
        }
    }
    catch(...)
    {
        delete [] edges;
        throw; //range test fails
    }
    auto * eStart = &edges[0];

    //2. Remove duplicate vertices, and (when closed) collinear edges ...
    auto * e = eStart;
    auto * eLoopStop = eStart;
    for(;;) {
        //nb: allows matching start and end points when not Closed ...
        if(e->curr == e->next->curr && (closed || e->next != eStart)) {
            if(e == e->next) break;
            if(e == eStart) eStart = e->next;
            e = RemoveEdge(e);
            eLoopStop = e;
            continue;
        }
        if(e->prev == e->next) break; //only two vertices
        else if(closed && SlopesEqual(e->prev->curr, e->curr, e->next->curr, m_useFullRange) &&
                (!m_preserveCollinear || !Pt2IsBetweenPt1AndPt3(e->prev->curr, e->curr, e->next->curr))) {
            //Collinear edges are allowed for open paths but in closed paths
            //the default is to merge adjacent collinear edges into a single edge.
            //However, if the PreserveCollinear property is enabled, only overlapping
            //collinear edges (ie spikes) will be removed from closed paths.
            if (e == eStart) eStart = e->next;
            e = RemoveEdge(e);
            e = e->prev;
            eLoopStop = e;
            continue;
        }
        e = e->next;
        if((e == eLoopStop) || (!closed && e->next == eStart)) break;
    }

    if((!closed && (e == e->next)) || (closed && (e->prev == e->next))) {
        delete [] edges;
        return false;
    }

    if(!closed) {
        m_hasOpenPaths = true;
        eStart->prev->outIdx = constant::skip;
    }

    //3. Do second stage of edge initialization ...
    e = eStart;
    do {
        InitEdge2(*e, polyTyp);
        e = e->next;
        if(isFlat && e->curr[1] != eStart->curr[1]) isFlat = false;
    } while(e != eStart);

    //4. Finally, add edge bounds to LocalMinima list ...

    //Totally flat paths must be handled differently when adding them
    //to LocalMinima list to avoid endless loops etc ...
    if(isFlat) {
        if(closed) {
            delete [] edges;
            return false;
        }
        e->prev->outIdx = constant::skip;
        typename MinimaList::value_type locMin;
        locMin.y = e->bot[1];
        locMin.boundL = nullptr;
        locMin.boundR = e;
        locMin.boundR->side = EdgeSide::Right;
        locMin.boundR->windDelta = 0;
        for (;;) {
            if(e->bot[0] != e->prev->top[0]) ReverseHorizontal(*e);
            if(e->next->outIdx == constant::skip) break;
            e->nextInLML = e->next;
            e = e->next;
        }
        m_minimaList.push_back(locMin);
        m_edges.push_back(edges);
        return true;
    }

    m_edges.push_back(edges);
    bool leftBoundIsForward;
    TEdge<num_type> * eMin = nullptr;

    //workaround to avoid an endless loop in the while loop below when
    //open paths have matching start and end points ...
    if(e->prev->bot == e->prev->top) e = e->next;

    for(;;) {
        e = FindNextLocMin(e);
        if (e == eMin) break;
        else if (!eMin) eMin = e;

        //e and e.Prev now share a local minima (left aligned if horizontal).
        //Compare their slopes to find which starts which bound ...
        typename MinimaList::value_type locMin;
        locMin.y = e->bot[1];
        if(math::LT(e->dx, e->prev->dx)) {
            locMin.boundL = e->prev;
            locMin.boundR = e;
            leftBoundIsForward = false; //Q.nextInLML = Q.prev
        }
        else {
            locMin.boundL = e;
            locMin.boundR = e->prev;
            leftBoundIsForward = true; //Q.nextInLML = Q.next
        }

        if(!closed) locMin.boundL->windDelta = 0;
        else if(locMin.boundL->next == locMin.boundR)
            locMin.boundL->windDelta = -1;
        else locMin.boundL->windDelta = 1;
        locMin.boundR->windDelta = -locMin.boundL->windDelta;

        e = ProcessBound(locMin.boundL, leftBoundIsForward);
        if (e->outIdx == constant::skip) e = ProcessBound(e, leftBoundIsForward);

        auto * e2 = ProcessBound(locMin.boundR, !leftBoundIsForward);
        if(e2->outIdx == constant::skip)
            e2 = ProcessBound(e2, !leftBoundIsForward);

        if(locMin.boundL->outIdx == constant::skip)
            locMin.boundL = nullptr;
        else if(locMin.boundR->outIdx == constant::skip)
            locMin.boundR = nullptr;
        m_minimaList.push_back(locMin);
        if (!leftBoundIsForward) e = e2;
    }
    return true;
}

template <typename num_type>
inline bool ClipperBase<num_type>::AddPaths(const Paths<num_type> & paths, PolyType polyTyp, bool closed)
{
    bool result = true;
    for(const auto & path : paths)
        result = result && AddPath(path, polyTyp, closed);
    return result;
}

template <typename num_type>
inline TEdge<num_type> * ClipperBase<num_type>::ProcessBound(TEdge<num_type> * e, bool nextIsForward)
{
    TEdge<num_type> * result = e;
    TEdge<num_type> * horz = nullptr;

    if(e->outIdx == constant::skip) {
        //if edges still remain in the current bound beyond the skip edge then
        //create another LocMin and call ProcessBound once more
        if(nextIsForward) {
            while (e->top[1] == e->next->bot[1]) e = e->next;
            //don't include top horizontals when parsing a bound a second time,
            //they will be contained in the opposite bound ...
            while (e != result && isHorizontal(*e)) e = e->prev;
        }
        else{
            while (e->top[1] == e->prev->bot[1]) e = e->prev;
            while (e != result && isHorizontal(*e)) e = e->next;
        }

        if(e == result) {
            if(nextIsForward) result = e->next;
            else result = e->prev;
        }
        else {
            //there are more edges in the bound beyond result starting with e
            if(nextIsForward)
                e = result->next;
            else
                e = result->prev;
            typename MinimaList::value_type locMin;
            locMin.y = e->bot[1];
            locMin.boundL = nullptr;
            locMin.boundR = e;
            e->windDelta = 0;
            result = ProcessBound(e, nextIsForward);
            m_minimaList.push_back(locMin);
        }
        return result;
    }

    TEdge<num_type> * eStart = nullptr;

    if(isHorizontal(*e)) {
        //We need to be careful with open paths because this may not be a
        //true local minima (ie e may be following a skip edge).
        //Also, consecutive horz. edges may start heading left before going right.
        if(nextIsForward)
            eStart = e->prev;
        else
            eStart = e->next;
        if(isHorizontal(*eStart)) {//ie an adjoining horizontal skip edge
            if(eStart->bot[0] != e->bot[0] && eStart->top[0] != e->bot[0])
                ReverseHorizontal(*e);
        }
        else if(eStart->bot[0] != e->bot[0])
            ReverseHorizontal(*e);
    }

     eStart = e;
    if(nextIsForward)
    {
        while (result->top[1] == result->next->bot[1] && result->next->outIdx != constant::skip)
            result = result->next;
        if(isHorizontal(*result) && result->next->outIdx != constant::skip) {
            //nb: at the top of a bound, horizontals are added to the bound
            //only when the preceding edge attaches to the horizontal's left vertex
            //unless a Skip edge is encountered when that becomes the top divide
            horz = result;
            while(isHorizontal(*horz->prev)) horz = horz->prev;
            if(horz->prev->top[0] > result->next->top[0]) result = horz->prev;
        }
        while(e != result) {
            e->nextInLML = e->next;
            if(isHorizontal(*e) && e != eStart &&
                e->bot[0] != e->prev->top[0]) ReverseHorizontal(*e);
            e = e->next;
        }
        if(isHorizontal(*e) && e != eStart && e->bot[0] != e->prev->top[0])
            ReverseHorizontal(*e);
        result = result->next; //move to the edge just beyond current bound
    }
    else {
        while(result->top[1] == result->prev->bot[1] && result->prev->outIdx != constant::skip)
            result = result->prev;
        if(isHorizontal(*result) && result->prev->outIdx != constant::skip) {
            horz = result;
            while(isHorizontal(*horz->next)) horz = horz->next;
            if(horz->next->top[0] == result->prev->top[0] ||
                horz->next->top[0] > result->prev->top[0])
                result = horz->next;
        }

        while(e != result) {
            e->nextInLML = e->prev;
            if(isHorizontal(*e) && e != eStart && e->bot[0] != e->next->top[0])
                ReverseHorizontal(*e);
            e = e->prev;
        }
        if(isHorizontal(*e) && e != eStart && e->bot[0] != e->next->top[0])
            ReverseHorizontal(*e);
        result = result->prev; //move to the edge just beyond current bound
    }

    return result;
}

template <typename num_type>
inline void ClipperBase<num_type>::InsertScanbeam(num_type y)
{
    m_scanbeam.push(y);
}

template <typename num_type>
inline bool ClipperBase<num_type>::PopScanbeam(num_type & y)
{
    if(m_scanbeam.empty()) return false;
    y = m_scanbeam.top();
    m_scanbeam.pop();
    while(!m_scanbeam.empty() && math::EQ(y, m_scanbeam.top())) { m_scanbeam.pop(); } // Pop duplicates.
    return true;  
}

template <typename num_type>
inline bool ClipperBase<num_type>::LocalMinimaPending()
{
    return m_currentLM != m_minimaList.end();
}

template <typename num_type>
inline bool ClipperBase<num_type>::PopLocalMinima(num_type y, const LocalMinimum<num_type> * & locMin)
{
    if(m_currentLM == m_minimaList.end() || math::NE(m_currentLM->y, y)) return false;
    locMin = &(*m_currentLM);
    ++m_currentLM;
    return true;
}

template <typename num_type>
inline OutRec<num_type> * ClipperBase<num_type>::CreateOutRec()
{
    m_polyOuts.push_back(new OutRec<num_type>);
    m_polyOuts.back()->idx = static_cast<int>(m_polyOuts.size()) - 1;
    return m_polyOuts.back();   
}

template <typename num_type>
inline void ClipperBase<num_type>::DisposeAllOutRecs()
{
    for(size_t i = 0; i < m_polyOuts.size(); ++i)
        DisposeOutRec(i);
    m_polyOuts.clear();
}


template <typename num_type>
inline void ClipperBase<num_type>::DisposeOutRec(size_t index)
{
    auto * outRec = m_polyOuts[index];
    if (outRec->pts) DisposeOutPts(outRec->pts);
    delete outRec;
    m_polyOuts[index] = nullptr;
}

template <typename num_type>
inline void ClipperBase<num_type>::SwapPositionsInAEL(TEdge<num_type> * edge1, TEdge<num_type> * edge2)
{
    //check that one or other edge hasn't already been removed from AEL ...
    if( edge1->nextInAEL == edge1->prevInAEL ||
        edge2->nextInAEL == edge2->prevInAEL ) return;

    if(edge1->nextInAEL == edge2) {
        auto * next = edge2->nextInAEL;
        if(next) next->prevInAEL = edge1;
        auto * prev = edge1->prevInAEL;
        if(prev) prev->nextInAEL = edge2;
        edge2->prevInAEL = prev;
        edge2->nextInAEL = edge1;
        edge1->prevInAEL = edge2;
        edge1->nextInAEL = next;
    }
    else if(edge2->nextInAEL == edge1) {
        auto * next = edge1->nextInAEL;
        if(next) next->prevInAEL = edge2;
        auto * prev = edge2->prevInAEL;
        if(prev) prev->nextInAEL = edge1;
        edge1->prevInAEL = prev;
        edge1->nextInAEL = edge2;
        edge2->prevInAEL = edge1;
        edge2->nextInAEL = next;
    }
    else {
        auto * next = edge1->nextInAEL;
        auto * prev = edge1->prevInAEL;
        edge1->nextInAEL = edge2->nextInAEL;
        if(edge1->nextInAEL) edge1->nextInAEL->prevInAEL = edge1;
        edge1->prevInAEL = edge2->prevInAEL;
        if(edge1->prevInAEL) edge1->prevInAEL->nextInAEL = edge1;
        edge2->nextInAEL = next;
        if(edge2->nextInAEL) edge2->nextInAEL->prevInAEL = edge2;
        edge2->prevInAEL = prev;
        if (edge2->prevInAEL) edge2->prevInAEL->nextInAEL = edge2;
    }

    if(!edge1->prevInAEL) m_activeEdges = edge1;
    else if(!edge2->prevInAEL) m_activeEdges = edge2;
}

template <typename num_type>
inline void ClipperBase<num_type>::DeleteFromAEL(TEdge<num_type> * e)
{
    auto * aelPrev = e->prevInAEL;
    auto * aelNext = e->nextInAEL;
    if (!aelPrev && !aelNext && e != m_activeEdges) return; //already deleted
    if (aelPrev) aelPrev->nextInAEL = aelNext;
    else m_activeEdges = aelNext;
    if (aelNext) aelNext->prevInAEL = aelPrev;
    e->nextInAEL = nullptr;
    e->prevInAEL = nullptr;
}

template <typename num_type>
inline void ClipperBase<num_type>::UpdateEdgeIntoAEL(TEdge<num_type> * & e)
{
    if(!e->nextInLML)
        throw Exception("UpdateEdgeIntoAEL: invalid call");

    e->nextInLML->outIdx = e->outIdx;
    auto * aelPrev = e->prevInAEL;
    auto * aelNext = e->nextInAEL;
    if(aelPrev) aelPrev->nextInAEL = e->nextInLML;
    else m_activeEdges = e->nextInLML;
    if(aelNext) aelNext->prevInAEL = e->nextInLML;
    e->nextInLML->side = e->side;
    e->nextInLML->windDelta = e->windDelta;
    e->nextInLML->windCnt = e->windCnt;
    e->nextInLML->windCnt2 = e->windCnt2;
    e = e->nextInLML;
    e->curr = e->bot;
    e->prevInAEL = aelPrev;
    e->nextInAEL = aelNext;
    if(!isHorizontal(*e)) InsertScanbeam(e->top[1]);
}

template <typename num_type>
class Clipper : public ClipperBase<num_type>
{
    using MaximaList = std::list<num_type>;
    using ClipperBase<num_type>::m_polyOuts;
    using ClipperBase<num_type>::m_activeEdges;
    using ClipperBase<num_type>::m_useFullRange;
    using ClipperBase<num_type>::m_preserveCollinear;
public:
    Clipper(int initOptions = 0);
    ~Clipper() = default;

    bool Execute(ClipType clipType, Paths<num_type> & solution, PolyFillType fillType);
    bool Execute(ClipType clipType, Paths<num_type> & solution, PolyFillType subjFillType, PolyFillType clipFillType);
    bool Execute(ClipType clipType, PolyTree<num_type> & polyTree, PolyFillType fillType = PolyFillType::EvenOdd);
    bool Execute(ClipType clipType, PolyTree<num_type> & polyTree, PolyFillType subjFillType, PolyFillType clipFillType);

private:
    bool ExecuteInternal();
    void SetWindingCount(TEdge<num_type> & edge);
    bool isEvenOddFillType(const TEdge<num_type> & edge) const;
    bool isEvenOddAltFillType(const TEdge<num_type> & edge) const;
    void InsertLocalMinimaIntoAEL(num_type botY);
    void InsertEdgeIntoAEL(TEdge<num_type> * edge, TEdge<num_type> * startEdge);
    void AddEdgeToSEL(TEdge<num_type> * edge);
    bool PopEdgeFromSEL(TEdge<num_type> * & edge);
    void CopyAELToSEL();
    void DeleteFromSEL(TEdge<num_type> * e);
    void SwapPositionsInSEL(TEdge<num_type> * edge1, TEdge<num_type> * edge2);
    bool isContributing(const TEdge<num_type> & edge) const;
    void DoMaxima(TEdge<num_type> * e);
    void ProcessHorizontals();
    void ProcessHorizontal(TEdge<num_type> * horzEdge);
    void AddLocalMaxPoly(TEdge<num_type> * e1, TEdge<num_type> * e2, const Point<num_type> & pt);
    OutPt<num_type> * AddLocalMinPoly(TEdge<num_type> * e1, TEdge<num_type> * e2, const Point<num_type> & pt);
    OutRec<num_type> * GetOutRec(size_t idx);
    void AppendPolygon(TEdge<num_type> * e1, TEdge<num_type> * e2);
    void IntersectEdges(TEdge<num_type> * e1, TEdge<num_type> * e2, Point<num_type> & pt);
    OutPt<num_type> * AddOutPt(TEdge<num_type> * e, const Point<num_type> & pt);
    OutPt<num_type> * GetLastOutPt(TEdge<num_type> * e);
    bool ProcessIntersections(const num_type topY);
    void BuildIntersectList(const num_type topY);
    void ProcessEdgesAtTopOfScanbeam(const num_type topY);
    void ProcessIntersectList();
    void BuildResult(Paths<num_type> & polys);
    void BuildResult2(PolyTree<num_type> & polytree);
    void SetHoleState(TEdge<num_type> * e, OutRec<num_type> * outrec);
    void DisposeIntersectNodes();
    void FixupOutPolygon(OutRec<num_type> & outrec);
    void FixupOutPolyline(OutRec<num_type> & outrec);
    bool FixupIntersectionOrder();
    void FixHoleLinkage(OutRec<num_type> & outrec);
    void AddJoin(OutPt<num_type> * op1, OutPt<num_type> * op2, Point<num_type> offPt);
    void AddGhostJoin(OutPt<num_type> * op, Point<num_type> offPt);
    void ClearJoins();
    void ClearGhostJoins();
    bool JoinPoints(Join<num_type> * j, OutRec<num_type> * outRec1, OutRec<num_type> * outRec2);
    void JoinCommonEdges();
    void DoSimplePolygons();
    void FixupFirstLefts1(OutRec<num_type> * oldOutRec, OutRec<num_type> * newOutRec);
    void FixupFirstLefts2(OutRec<num_type> * innerOutRec, OutRec<num_type> * outerOutRec);
    void FixupFirstLefts3(OutRec<num_type> * oldOutRec, OutRec<num_type> * newOutRec);

private:
    JoinList<num_type> m_joins;
    JoinList<num_type> m_ghostJoins;
    IntersectList<num_type> m_intersectList;
    ClipType m_clipType;
    MaximaList m_maxima;
    TEdge<num_type> * m_sortedEdges;
    PolyFillType m_clipFillType;
    PolyFillType m_subjFillType;
    bool m_executeLocked;
    bool m_reverseOutput;
    bool m_usingPolyTree;
    bool m_strictSimple;
};

template <typename num_type>
inline Clipper<num_type>::Clipper(int initOptions) : ClipperBase<num_type>()
{
    m_executeLocked = false;
    ClipperBase<num_type>::m_useFullRange = false;
    m_reverseOutput = (initOptions & static_cast<int>(InitOptions::ReverseSolution)) != 0;
    m_strictSimple = (initOptions & static_cast<int>(InitOptions::StrictlySimple)) != 0;
    ClipperBase<num_type>::m_preserveCollinear = (initOptions & static_cast<int>(InitOptions::PreserveCollinear)) != 0;
    ClipperBase<num_type>::m_hasOpenPaths = false;
}

template <typename num_type>
inline bool Clipper<num_type>::Execute(ClipType clipType, Paths<num_type> & solution, PolyFillType fillType)
{
    return Execute(clipType, solution, fillType, fillType);
}

template <typename num_type>
inline bool Clipper<num_type>::Execute(ClipType clipType, Paths<num_type> & solution, PolyFillType subjFillType, PolyFillType clipFillType)
{
    if(m_executeLocked) return false;
    if(ClipperBase<num_type>::m_hasOpenPaths)
        throw Exception("Error: PolyTree struct is needed for open path clipping.");
    m_executeLocked = true;
    solution.clear();
    m_subjFillType = subjFillType;
    m_clipFillType = clipFillType;
    m_clipType = clipType;
    m_usingPolyTree = false;
    auto res = ExecuteInternal();
    if(res) BuildResult(solution);
    ClipperBase<num_type>::DisposeAllOutRecs();
    m_executeLocked = false;
    return res;
}

template <typename num_type>
inline bool Clipper<num_type>::Execute(ClipType clipType, PolyTree<num_type> & polyTree, PolyFillType fillType)
{
    return Execute(clipType, polyTree, fillType, fillType);
}

template <typename num_type>
inline bool Clipper<num_type>::Execute(ClipType clipType, PolyTree<num_type> & polyTree, PolyFillType subjFillType, PolyFillType clipFillType)
{
    if(m_executeLocked) return false;
    m_executeLocked = true;
    m_subjFillType = subjFillType;
    m_clipFillType = clipFillType;
    m_clipType = clipType;
    m_usingPolyTree = true;
    auto res = ExecuteInternal();
    if (res) BuildResult2(polyTree);
    ClipperBase<num_type>::DisposeAllOutRecs();
    m_executeLocked = false;
    return res;
}

template <typename num_type>
inline bool Clipper<num_type>::ExecuteInternal()
{
    bool succeeded = true;
    try {
        ClipperBase<num_type>::Reset();
        m_maxima.clear();
        m_sortedEdges = nullptr;

        succeeded = true;
        num_type botY, topY;
        if(!ClipperBase<num_type>::PopScanbeam(botY)) return false;
        InsertLocalMinimaIntoAEL(botY);
        while(ClipperBase<num_type>::PopScanbeam(topY) || ClipperBase<num_type>::LocalMinimaPending()) {
            ProcessHorizontals();
            ClearGhostJoins();
            if(!ProcessIntersections(topY)) {
                succeeded = false;
                break;
            }
            ProcessEdgesAtTopOfScanbeam(topY);
            botY = topY;
            InsertLocalMinimaIntoAEL(botY);
        }
    }
    catch(...) {
        succeeded = false;
    }

    if (succeeded) {
        //fix orientations ...
        for(auto * outRec : m_polyOuts) {
            if (!outRec->pts || outRec->isOpen) continue;
            if ((outRec->isHole ^ m_reverseOutput) == math::isPositive(Area(*outRec)))
                ReversePolyPtLinks(outRec->pts);
        }

        if(!m_joins.empty()) JoinCommonEdges();

        //unfortunately FixupOutPolygon() must be done after JoinCommonEdges()
        for (auto * outRec : m_polyOuts) {
            if (!outRec->pts) continue;
            if (outRec->isOpen)
                FixupOutPolyline(*outRec);
            else FixupOutPolygon(*outRec);
        }

        if(m_strictSimple)
            DoSimplePolygons();
    }

    ClearJoins();
    ClearGhostJoins();
    return succeeded;
}

template <typename num_type>
inline void Clipper<num_type>::BuildIntersectList(const num_type topY)
{
    if (!m_activeEdges ) return;

    //prepare for sorting ...
    auto * e = m_activeEdges;
    m_sortedEdges = e;
    while(e) {
        e->prevInSEL = e->prevInAEL;
        e->nextInSEL = e->nextInAEL;
        e->curr[0] = TopX(*e, topY);
        e = e->nextInAEL;
    }

    //bubblesort ...
    bool isModified;
    do {
        isModified = false;
        e = m_sortedEdges;
        while(e->nextInSEL) {
            auto * eNext = e->nextInSEL;
            Point<num_type> pt;
            if(math::GT(e->curr[0], eNext->curr[0])) {
                IntersectPoint(*e, *eNext, pt);
                if(math::LT(pt[1], topY))
                    pt = Point<num_type>(TopX(*e, topY), topY);
                auto * newNode = new IntersectNode<num_type>;
                newNode->edge1 = e;
                newNode->edge2 = eNext;
                newNode->point = pt;
                m_intersectList.push_back(newNode);

                SwapPositionsInSEL(e, eNext);
                isModified = true;
            }
            else e = eNext;
        }
        if(e->prevInSEL)
            e->prevInSEL->nextInSEL = nullptr;
        else break;
    } while(isModified);
    m_sortedEdges = nullptr; //important    
}

template <typename num_type>
inline void Clipper<num_type>::ProcessEdgesAtTopOfScanbeam(const num_type topY)
{
    auto * e = m_activeEdges;
    while(e) {
        //1. process maxima, treating them as if they're 'bent' horizontal edges,
        //   but exclude maxima with horizontal edges. nb: e can't be a horizontal.
        bool isMaximaEdge = isMaxima(e, topY);

        if(isMaximaEdge) {
            auto * eMaxPair = GetMaximaPairEx(e);
            isMaximaEdge = !eMaxPair || !isHorizontal(*eMaxPair);
        }

        if(isMaximaEdge) {
            if (m_strictSimple)
                m_maxima.push_back(e->top[0]);
            auto * ePrev = e->prevInAEL;
            DoMaxima(e);
            if(!ePrev)
                e = m_activeEdges;
            else e = ePrev->nextInAEL;
        }
        else {
            //2. promote horizontal edges, otherwise update Curr[0] and Curr[1] ...
            if(isIntermediate(e, topY) && isHorizontal(*e->nextInLML)) {
                ClipperBase<num_type>::UpdateEdgeIntoAEL(e);
                if(e->outIdx >= 0) 
                    AddOutPt(e, e->bot);
                AddEdgeToSEL(e);
            }
            else {
                e->curr[0] = TopX(*e, topY);
                e->curr[1] = topY;
            }

            //When StrictlySimple and 'e' is being touched by another edge, then
            //make sure both edges have a vertex here ...
            if(m_strictSimple) {
                auto * ePrev = e->prevInAEL;
                if(e->outIdx >= 0 && e->windDelta != 0 && ePrev && ePrev->outIdx >= 0 &&
                    ePrev->curr[0] == e->curr[0] && ePrev->windDelta != 0) {
                    auto pt = e->curr;
                    auto * op1 = AddOutPt(ePrev, pt);
                    auto * op2 = AddOutPt(e, pt);
                    AddJoin(op1, op2, pt); //StrictlySimple (type-3) join
                }
            }
            e = e->nextInAEL;
        }
    }

    //3. Process horizontals at the Top of the scanbeam ...
    m_maxima.sort();
    ProcessHorizontals();
    m_maxima.clear();

    //4. Promote intermediate vertices ...
    e = m_activeEdges;
    while(e) {
        if(isIntermediate(e, topY)) {
            OutPt<num_type> * op1 = nullptr;
            if(e->outIdx >= 0)
                op1 = AddOutPt(e, e->top);
            ClipperBase<num_type>::UpdateEdgeIntoAEL(e);

            //if output polygons share an edge, they'll need joining later ...
            auto * ePrev = e->prevInAEL;
            auto * eNext = e->nextInAEL;
            if (ePrev && ePrev->curr[0] == e->bot[0] &&
                ePrev->curr[1] == e->bot[1] && op1 &&
                ePrev->outIdx >= 0 && ePrev->curr[1] > ePrev->top[1] &&
                SlopesEqual(e->curr, e->top, ePrev->curr, ePrev->top, m_useFullRange) &&
                e->windDelta != 0 && ePrev->windDelta != 0) {
                auto * op2 = AddOutPt(ePrev, e->bot);
                AddJoin(op1, op2, e->top);
            }
            else if (eNext && eNext->curr[0] == e->bot[0] &&
                eNext->curr[1] == e->bot[1] && op1 &&
                eNext->outIdx >= 0 && eNext->curr[1] > eNext->top[1] &&
                SlopesEqual(e->curr, e->top, eNext->curr, eNext->top, m_useFullRange) &&
                e->windDelta != 0 && eNext->windDelta != 0) {
                auto * op2 = AddOutPt(eNext, e->bot);
                AddJoin(op1, op2, e->top);
            }
        }
        e = e->nextInAEL;
    }
}

template <typename num_type>
inline void Clipper<num_type>::ProcessIntersectList()
{
    for(auto * node : m_intersectList) {
        IntersectEdges(node->edge1, node->edge2, node->point);
        ClipperBase<num_type>::SwapPositionsInAEL(node->edge1 , node->edge2);
        delete node;
    }
    m_intersectList.clear();
}

template <typename num_type>
inline void Clipper<num_type>::BuildResult(Paths<num_type> & polys)
{
    polys.reserve(m_polyOuts.size());
    for (auto * outRec : m_polyOuts) {
        if(!outRec->pts) continue;

        Path<num_type> path;
        auto * p = outRec->pts->prev;
        auto size = PointCount(p);
        if (size < 2) continue;
        path.reserve(size);
        for(auto i = 0; i < size; ++i) {
            path.push_back(p->pt);
            p = p->prev;
        }
        polys.emplace_back(std::move(path));
    }
}

template <typename num_type>
inline void Clipper<num_type>::BuildResult2(PolyTree<num_type> & polytree)
{
    polytree.Clear();
    polytree.m_allNodes.reserve(m_polyOuts.size());
    //add each output polygon/contour to polytree ...
    for (auto * outRec : m_polyOuts) {
        auto size = PointCount(outRec->pts);
        if((outRec->isOpen && size < 2) || (!outRec->isOpen && size < 3)) continue;
        FixHoleLinkage(*outRec);
        auto * pn = new PolyNode<num_type>;
        //nb: polytree takes ownership of all the PolyNodes
        polytree.m_allNodes.push_back(pn);
        outRec->polyNd = pn;
        pn->parent = nullptr;
        pn->m_index = 0;
        pn->contour.reserve(size);
        auto * op = outRec->pts->prev;
        for (auto i = 0; i < size; ++i) {
            pn->contour.push_back(op->pt);
            op = op->prev;
        }
    }

    //fixup PolyNode links etc ...
    polytree.children.reserve(m_polyOuts.size());
    for (auto * outRec : m_polyOuts) {
        if(!outRec->polyNd) continue;
        if(outRec->isOpen) {
            outRec->polyNd->m_isOpen = true;
            polytree.AddChild(*outRec->polyNd);
        }
        else if(outRec->firstLeft && outRec->firstLeft->polyNd)
            outRec->firstLeft->polyNd->AddChild(*outRec->polyNd);
        else polytree.AddChild(*outRec->polyNd);
    }
}

template <typename num_type>
inline void Clipper<num_type>::SetWindingCount(TEdge<num_type> & edge)
{
    auto * e = edge.prevInAEL;
    //find the edge of the same polytype that immediately preceeds 'edge' in AEL
    while(e && ((e->polyTyp != edge.polyTyp) || (e->windDelta == 0))) e = e->prevInAEL;
    if(!e) {
        if (edge.windDelta == 0) {
            PolyFillType pft = edge.polyTyp == PolyType::Subject ? m_subjFillType : m_clipFillType;
            edge.windCnt = pft == PolyFillType::Negative ? -1 : 1;
        }
        else edge.windCnt = edge.windDelta;
        edge.windCnt2 = 0;
        e = m_activeEdges; //ie get ready to calc WindCnt2
    }
    else if (edge.windDelta == 0 && m_clipType != ClipType::Union) {
        edge.windCnt = 1;
        edge.windCnt2 = e->windCnt2;
        e = e->nextInAEL; //ie get ready to calc WindCnt2
    }
    else if(isEvenOddFillType(edge))
    {
        //EvenOdd filling ...
        if(edge.windDelta == 0) {
            //are we inside a subj polygon ...
            bool inside = true;
            auto * e2 = e->prevInAEL;
            while(e2) {
                if(e2->polyTyp == e->polyTyp && e2->windDelta != 0)
                    inside = !inside;
                e2 = e2->prevInAEL;
            }
            edge.windCnt = inside ? 0 : 1;
        }
        else edge.windCnt = edge.windDelta;
        edge.windCnt2 = e->windCnt2;
        e = e->nextInAEL; //ie get ready to calc WindCnt2
    }
    else
    {
        //nonZero, Positive or Negative filling ...
        if(e->windCnt * e->windDelta < 0) {
            //prev edge is 'decreasing' WindCount (WC) toward zero
            //so we're outside the previous polygon ...
            if (std::abs(e->windCnt) > 1) {
                //outside prev poly but still inside another.
                //when reversing direction of prev poly use the same WC
                if (e->windDelta * edge.windDelta < 0) edge.windCnt = e->windCnt;
                //otherwise continue to 'decrease' WC ...
                else edge.windCnt = e->windCnt + edge.windDelta;
            }
            else
                //now outside all polys of same polytype so set own WC ...
                edge.windCnt = edge.windDelta == 0 ? 1 : edge.windDelta;
        }
        else {
            //prev edge is 'increasing' WindCount (WC) away from zero
            //so we're inside the previous polygon ...
            if(edge.windDelta == 0)
                edge.windCnt = e->windCnt < 0 ? e->windCnt - 1 : e->windCnt + 1;
            //if wind direction is reversing prev then use same WC
            else if(e->windDelta * edge.windDelta < 0) edge.windCnt = e->windCnt;
            //otherwise add to WC ...
            else edge.windCnt = e->windCnt + edge.windDelta;
        }
        edge.windCnt2 = e->windCnt2;
        e = e->nextInAEL; //ie get ready to calc WindCnt2
    }

    //update WindCnt2 ...
    if(isEvenOddAltFillType(edge)) {
        //EvenOdd filling ...
        while(e != &edge) {
            if(e->windDelta != 0)
                edge.windCnt2 = edge.windCnt2 == 0 ? 1 : 0;
            e = e->nextInAEL;
        }
    }
    else {
        //nonZero, Positive or Negative filling ...
        while(e != &edge) {
            edge.windCnt2 += e->windDelta;
            e = e->nextInAEL;
        }
    }   
}

template <typename num_type>
inline bool Clipper<num_type>::isEvenOddFillType(const TEdge<num_type> & edge) const
{
    if (edge.polyTyp == PolyType::Subject)
        return m_subjFillType == PolyFillType::EvenOdd;
    else return m_clipFillType == PolyFillType::EvenOdd;
}

template <typename num_type>
inline bool Clipper<num_type>::isEvenOddAltFillType(const TEdge<num_type> & edge) const
{
    if (edge.polyTyp == PolyType::Subject)
        return m_clipFillType == PolyFillType::EvenOdd;
    else return m_subjFillType == PolyFillType::EvenOdd;
}

template <typename num_type>
inline void Clipper<num_type>::InsertLocalMinimaIntoAEL(num_type botY)
{
    const LocalMinimum<num_type> * lm = nullptr;
    while(ClipperBase<num_type>::PopLocalMinima(botY, lm)){
        auto * lb = lm->boundL;
        auto * rb = lm->boundR;

        OutPt<num_type> * op1 = nullptr;
        if(!lb) {
            //nb: don't insert LB into either AEL or SEL
            InsertEdgeIntoAEL(rb, nullptr);
            SetWindingCount(*rb);
            if(isContributing(*rb))
                op1 = AddOutPt(rb, rb->bot);
        }
        else if(!rb) {
            InsertEdgeIntoAEL(lb, nullptr);
            SetWindingCount(*lb);
            if (isContributing(*lb))
                op1 = AddOutPt(lb, lb->bot);
            ClipperBase<num_type>::InsertScanbeam(lb->top[1]);
        }
        else {
            InsertEdgeIntoAEL(lb, nullptr);
            InsertEdgeIntoAEL(rb, lb);
            SetWindingCount(*lb);
            rb->windCnt = lb->windCnt;
            rb->windCnt2 = lb->windCnt2;
            if(isContributing(*lb))
                op1 = AddLocalMinPoly(lb, rb, lb->bot);
            ClipperBase<num_type>::InsertScanbeam(lb->top[1]);
        }

        if(rb) {
            if(isHorizontal(*rb)) {
                AddEdgeToSEL(rb);
                if(rb->nextInLML)
                    ClipperBase<num_type>::InsertScanbeam(rb->nextInLML->top[1]);
            }
            else ClipperBase<num_type>::InsertScanbeam(rb->top[1]);
        }

        if(!lb || !rb) continue;

        //if any output polygons share an edge, they'll need joining later ...
        if(op1 && isHorizontal(*rb) && !m_ghostJoins.empty() && (rb->windDelta != 0)) {
            for(auto * jr : m_ghostJoins) {
                //if the horizontal Rb and a 'ghost' horizontal overlap, then convert
                //the 'ghost' join to a real join ready for later ...
                if(HorzSegmentsOverlap(jr->outPt1->pt[0], jr->offPt[0], rb->bot[0], rb->top[0]))
                    AddJoin(jr->outPt1, op1, jr->offPt);
            }
        }

        if(lb->outIdx >= 0 && lb->prevInAEL &&
            lb->prevInAEL->curr[0] == lb->bot[0] && lb->prevInAEL->outIdx >= 0 &&
            SlopesEqual(lb->prevInAEL->bot, lb->prevInAEL->top, lb->curr, lb->top, m_useFullRange) &&
            lb->windDelta != 0 && lb->prevInAEL->windDelta != 0) {
            auto * op2 = AddOutPt(lb->prevInAEL, lb->bot);
            AddJoin(op1, op2, lb->top);
        }

        if(lb->nextInAEL != rb) {
            if(rb->outIdx >= 0 && rb->prevInAEL->outIdx >= 0 &&
                SlopesEqual(rb->prevInAEL->curr, rb->prevInAEL->top, rb->curr, rb->top, m_useFullRange) &&
                rb->windDelta != 0 && rb->prevInAEL->windDelta != 0) {
                auto * op2 = AddOutPt(rb->prevInAEL, rb->bot);
                AddJoin(op1, op2, rb->top);
            }

            auto * e = lb->nextInAEL;
            if(e) {
                while(e != rb) {
                    //nb: For calculating winding counts etc, IntersectEdges() assumes
                    //that param1 will be to the Right of param2 ABOVE the intersection ...
                    IntersectEdges(rb ,e ,lb->curr); //order important here
                    e = e->nextInAEL;
                }
            }
        }
    }
}

template <typename num_type>
inline void Clipper<num_type>::SwapPositionsInSEL(TEdge<num_type> * edge1, TEdge<num_type> * edge2)
{
    if(!(edge1->nextInSEL) && !(edge1->prevInSEL)) return;
    if(!(edge2->nextInSEL) && !(edge2->prevInSEL)) return;

    if(edge1->nextInSEL == edge2 ) {
        auto * next = edge2->nextInSEL;
        if(next) next->prevInSEL = edge1;
        auto * prev = edge1->prevInSEL;
        if(prev) prev->nextInSEL = edge2;
        edge2->prevInSEL = prev;
        edge2->nextInSEL = edge1;
        edge1->prevInSEL = edge2;
        edge1->nextInSEL = next;
    }
    else if(edge2->nextInSEL == edge1) {
        auto * next = edge1->nextInSEL;
        if(next) next->prevInSEL = edge2;
        auto * prev = edge2->prevInSEL;
        if(prev) prev->nextInSEL = edge1;
        edge1->prevInSEL = prev;
        edge1->nextInSEL = edge2;
        edge2->prevInSEL = edge1;
        edge2->nextInSEL = next;
    }
    else {
        auto * next = edge1->nextInSEL;
        auto * prev = edge1->prevInSEL;
        edge1->nextInSEL = edge2->nextInSEL;
        if(edge1->nextInSEL) edge1->nextInSEL->prevInSEL = edge1;
        edge1->prevInSEL = edge2->prevInSEL;
        if(edge1->prevInSEL) edge1->prevInSEL->nextInSEL = edge1;
        edge2->nextInSEL = next;
        if(edge2->nextInSEL) edge2->nextInSEL->prevInSEL = edge2;
        edge2->prevInSEL = prev;
        if(edge2->prevInSEL) edge2->prevInSEL->nextInSEL = edge2;
    }

    if(!edge1->prevInSEL) m_sortedEdges = edge1;
    else if(!edge2->prevInSEL) m_sortedEdges = edge2;
}

template <typename num_type>
inline bool Clipper<num_type>::isContributing(const TEdge<num_type> & edge) const
{
    PolyFillType pft1, pft2;
    if(edge.polyTyp == PolyType::Subject) {
        pft1 = m_subjFillType;
        pft2 = m_clipFillType;
    }
    else {
        pft1 = m_clipFillType;
        pft2 = m_subjFillType;
    }

    switch(pft1) {
        case PolyFillType::EvenOdd :
            //return false if a subj line has been flagged as inside a subj polygon
            if(edge.windDelta == 0 && edge.windCnt != 1) return false;
            break;
        case PolyFillType::NonZero :
            if(std::abs(edge.windCnt) != 1) return false;
            break;
        case PolyFillType::Positive :
            if(edge.windCnt != 1) return false;
            break;
        default : //Negative
            if(edge.windCnt != -1) return false;
    }

    switch(m_clipType) {
        case ClipType::Intersection :
            switch(pft2)
            {
                case PolyFillType::EvenOdd :
                case PolyFillType::NonZero :
                    return edge.windCnt2 != 0;
                case PolyFillType::Positive :
                    return edge.windCnt2 > 0;
                default:
                    return edge.windCnt2 < 0;
            }
            break;
        case ClipType::Union :
            switch(pft2) {
                case PolyFillType::EvenOdd :
                case PolyFillType::NonZero :
                    return edge.windCnt2 == 0;
                case PolyFillType::Positive :
                    return edge.windCnt2 <= 0;
                default :
                    return edge.windCnt2 >= 0;
            }
            break;
        case ClipType::Difference :
            if (edge.polyTyp == PolyType::Subject) {
                switch(pft2) {
                    case PolyFillType::EvenOdd :
                    case PolyFillType::NonZero :
                        return edge.windCnt2 == 0;
                    case PolyFillType::Positive :
                        return edge.windCnt2 <= 0;
                    default :
                        return edge.windCnt2 >= 0;
                }
            }
            else {
                switch(pft2) {
                    case PolyFillType::EvenOdd :
                    case PolyFillType::NonZero :
                        return edge.windCnt2 != 0;
                    case PolyFillType::Positive :
                        return edge.windCnt2 > 0;
                    default :
                        return edge.windCnt2 < 0;
                }
            }
            break;
        case ClipType::Xor :
            if(edge.windDelta == 0) {//XOr always contributing unless open
                switch(pft2) {
                    case PolyFillType::EvenOdd :
                    case PolyFillType::NonZero :
                        return edge.windCnt2 == 0;
                    case PolyFillType::Positive :
                        return edge.windCnt2 <= 0;
                    default:
                        return edge.windCnt2 >= 0;
                }
            }
            else return true;
            break;
        default:
            return true;
    }
}

template <typename num_type>
inline void Clipper<num_type>::DoMaxima(TEdge<num_type> * e)
{
    auto * eMaxPair = GetMaximaPairEx(e);
    if(!eMaxPair) {
        if(e->outIdx >= 0) AddOutPt(e, e->top);
        ClipperBase<num_type>::DeleteFromAEL(e);
        return;
    }

    auto * eNext = e->nextInAEL;
    while(eNext && eNext != eMaxPair) {
        IntersectEdges(e, eNext, e->top);
        ClipperBase<num_type>::SwapPositionsInAEL(e, eNext);
        eNext = e->nextInAEL;
    }

    if(e->outIdx == constant::unassigned && eMaxPair->outIdx == constant::unassigned) {
        ClipperBase<num_type>::DeleteFromAEL(e);
        ClipperBase<num_type>::DeleteFromAEL(eMaxPair);
    }
    else if(e->outIdx >= 0 && eMaxPair->outIdx >= 0) {
        if (e->outIdx >= 0)
            AddLocalMaxPoly(e, eMaxPair, e->top);
        ClipperBase<num_type>::DeleteFromAEL(e);
        ClipperBase<num_type>::DeleteFromAEL(eMaxPair);
    }
    else if(e->windDelta == 0) {
        if(e->outIdx >= 0) {
            AddOutPt(e, e->top);
            e->outIdx = constant::unassigned;
        }
        ClipperBase<num_type>::DeleteFromAEL(e);

        if(eMaxPair->outIdx >= 0) {
            AddOutPt(eMaxPair, e->top);
            eMaxPair->outIdx = constant::unassigned;
        }
        ClipperBase<num_type>::DeleteFromAEL(eMaxPair);
    }
    else throw Exception("DoMaxima error"); 
}

template <typename num_type>
inline void Clipper<num_type>::ProcessHorizontals()
{
    TEdge<num_type> * horzEdge;
    while(PopEdgeFromSEL(horzEdge))
        ProcessHorizontal(horzEdge);
}

template <typename num_type>
inline void Clipper<num_type>::ProcessHorizontal(TEdge<num_type> * horzEdge)
{
    Direction dir;
    num_type horzLeft, horzRight;
    bool isOpen = horzEdge->windDelta == 0;

    GetHorzDirection(*horzEdge, dir, horzLeft, horzRight);

    TEdge<num_type> * eLastHorz = horzEdge, * eMaxPair = nullptr;
    while(eLastHorz->nextInLML && isHorizontal(*eLastHorz->nextInLML))
        eLastHorz = eLastHorz->nextInLML;
    if(!eLastHorz->nextInLML)
        eMaxPair = GetMaximaPair(eLastHorz);

    typename MaximaList::const_iterator maxIt;
    typename MaximaList::const_reverse_iterator maxRit;
    if(m_maxima.size() > 0) {
        //get the first maxima in range (X) ...
        if(dir == Direction::LeftToRight) {
            maxIt = m_maxima.begin();
            while(maxIt != m_maxima.end() && math::LE(*maxIt, horzEdge->bot[0])) maxIt++;
            if(maxIt != m_maxima.end() && math::GE(*maxIt, eLastHorz->top[0])) maxIt = m_maxima.end();
        }
        else {
            maxRit = m_maxima.rbegin();
            while(maxRit != m_maxima.rend() && math::GT(*maxRit, horzEdge->bot[0])) maxRit++;
            if(maxRit != m_maxima.rend() && math::LE(*maxRit, eLastHorz->top[0])) maxRit = m_maxima.rend();
        }
    }

    OutPt<num_type> * op1 = nullptr;

    for(;;) {//loop through consec. horizontal edges
        bool isLastHorz = horzEdge == eLastHorz;
        auto * e = GetNextInAEL(horzEdge, dir);
        while(e) {
            //this code block inserts extra coords into horizontal edges (in output
            //polygons) whereever maxima touch these horizontal edges. This helps
            //'simplifying' polygons (ie if the Simplify property is set).
            if(m_maxima.size() > 0) {
                if(dir == Direction::LeftToRight) {
                    while(maxIt != m_maxima.end() && math::LT(*maxIt, e->curr[0])) {
                        if(horzEdge->outIdx >= 0 && !isOpen)
                            AddOutPt(horzEdge, Point<num_type>(*maxIt, horzEdge->bot[1]));
                        maxIt++;
                    }
                }
                else {
                    while(maxRit != m_maxima.rend() && math::GT(*maxRit, e->curr[0])) {
                        if(horzEdge->outIdx >= 0 && !isOpen)
                            AddOutPt(horzEdge, Point<num_type>(*maxRit, horzEdge->bot[1]));
                        maxRit++;
                    }
                }
            }

            if((dir == Direction::LeftToRight && math::GT(e->curr[0], horzRight)) ||
                (dir == Direction::RightToLeft && math::LT(e->curr[0], horzLeft))) break;

            //Also break if we've got to the end of an intermediate horizontal edge ...
            //nb: Smaller Dx's are to the right of larger Dx's ABOVE the horizontal.
            if(math::EQ(e->curr[0], horzEdge->top[0]) && horzEdge->nextInLML &&
                math::LT(e->dx, horzEdge->nextInLML->dx)) break;

            if(horzEdge->outIdx >= 0 && !isOpen) {//note: may be done multiple times
                op1 = AddOutPt(horzEdge, e->curr);
                auto * eNextHorz = m_sortedEdges;
                while (eNextHorz) {
                    if (eNextHorz->outIdx >= 0 &&
                        HorzSegmentsOverlap(horzEdge->bot[0], horzEdge->top[0], eNextHorz->bot[0], eNextHorz->top[0])) {
                        auto * op2 = GetLastOutPt(eNextHorz);
                        AddJoin(op2, op1, eNextHorz->top);
                    }
                    eNextHorz = eNextHorz->nextInSEL;
                }
                AddGhostJoin(op1, horzEdge->bot);
            }

            //OK, so far we're still in range of the horizontal Edge  but make sure
            //we're at the last of consec. horizontals when matching with eMaxPair
            if(e == eMaxPair && isLastHorz) {
                if(horzEdge->outIdx >= 0)
                    AddLocalMaxPoly(horzEdge, eMaxPair, horzEdge->top);
                ClipperBase<num_type>::DeleteFromAEL(horzEdge);
                ClipperBase<num_type>::DeleteFromAEL(eMaxPair);
                return;
            }

            if(dir == Direction::LeftToRight) {
                auto pt = Point<num_type>(e->curr[0], horzEdge->curr[1]);
                IntersectEdges(horzEdge, e, pt);
            }
            else {
                auto pt = Point<num_type>(e->curr[0], horzEdge->curr[1]);
                IntersectEdges(e, horzEdge, pt);
            }
            auto * eNext = GetNextInAEL(e, dir);
            ClipperBase<num_type>::SwapPositionsInAEL(horzEdge, e);
            e = eNext;
        }

        //Break out of loop if HorzEdge.NextInLML is not also horizontal ...
        if(!horzEdge->nextInLML || !isHorizontal(*horzEdge->nextInLML)) break;

        ClipperBase<num_type>::UpdateEdgeIntoAEL(horzEdge);
        if(horzEdge->outIdx >= 0)
            AddOutPt(horzEdge, horzEdge->bot);
        GetHorzDirection(*horzEdge, dir, horzLeft, horzRight);

    }

    if(horzEdge->outIdx >= 0 && !op1) {
        op1 = GetLastOutPt(horzEdge);
        auto * eNextHorz = m_sortedEdges;
        while(eNextHorz) {
            if(eNextHorz->outIdx >= 0 &&
                HorzSegmentsOverlap(horzEdge->bot[0], horzEdge->top[0], eNextHorz->bot[0], eNextHorz->top[0])) {
                auto * op2 = GetLastOutPt(eNextHorz);
                AddJoin(op2, op1, eNextHorz->top);
            }
            eNextHorz = eNextHorz->nextInSEL;
        }
        AddGhostJoin(op1, horzEdge->top);
    }

    if (horzEdge->nextInLML) {
        if(horzEdge->outIdx >= 0) {
            op1 = AddOutPt(horzEdge, horzEdge->top);
            ClipperBase<num_type>::UpdateEdgeIntoAEL(horzEdge);
            if(horzEdge->windDelta == 0) return;
            //nb: HorzEdge is no longer horizontal here
            auto * ePrev = horzEdge->prevInAEL;
            auto * eNext = horzEdge->nextInAEL;
            if(ePrev &&
                math::EQ(ePrev->curr[0], horzEdge->bot[0]) && 
                math::EQ(ePrev->curr[1], horzEdge->bot[1]) &&
                ePrev->windDelta != 0 && ePrev->outIdx >= 0 &&
                math::GT(ePrev->curr[1], ePrev->top[1]) && SlopesEqual(*horzEdge, *ePrev, m_useFullRange)) {
                auto * op2 = AddOutPt(ePrev, horzEdge->bot);
                AddJoin(op1, op2, horzEdge->top);
            }
            else if(eNext &&
                math::EQ(eNext->curr[0], horzEdge->bot[0]) &&
                math::EQ(eNext->curr[1], horzEdge->bot[1]) &&
                eNext->windDelta != 0 && eNext->outIdx >= 0 &&
                math::GT(eNext->curr[1], eNext->top[1]) && SlopesEqual(*horzEdge, *eNext, m_useFullRange)) {
                auto * op2 = AddOutPt(eNext, horzEdge->bot);
                AddJoin(op1, op2, horzEdge->top);
            }
        }
        else ClipperBase<num_type>::UpdateEdgeIntoAEL(horzEdge);
    }
    else {
        if(horzEdge->outIdx >= 0)
            AddOutPt(horzEdge, horzEdge->top);
        ClipperBase<num_type>::DeleteFromAEL(horzEdge);
    } 
}

template <typename num_type>
inline void Clipper<num_type>::AddLocalMaxPoly(TEdge<num_type> * e1, TEdge<num_type> * e2, const Point<num_type> & pt)
{
    AddOutPt(e1, pt);
    if(e2->windDelta == 0) AddOutPt(e2, pt);
    if(e1->outIdx == e2->outIdx ) {
        e1->outIdx = constant::unassigned;
        e2->outIdx = constant::unassigned;
    }
    else if(e1->outIdx < e2->outIdx)
        AppendPolygon(e1, e2);
    else AppendPolygon(e2, e1);   
}

template <typename num_type>
inline OutPt<num_type> * Clipper<num_type>::AddLocalMinPoly(TEdge<num_type> * e1, TEdge<num_type> * e2, const Point<num_type> & pt)
{
    OutPt<num_type> * result = nullptr;
    TEdge<num_type> * e = nullptr, * prevE = nullptr;
    if(isHorizontal(*e2) || math::GT(e1->dx, e2->dx)) {
        result = AddOutPt(e1, pt);
        e2->outIdx = e1->outIdx;
        e1->side = EdgeSide::Left;
        e2->side = EdgeSide::Right;
        e = e1;
        if(e->prevInAEL == e2)
            prevE = e2->prevInAEL;
        else prevE = e->prevInAEL;
    }
    else {
        result = AddOutPt(e2, pt);
        e1->outIdx = e2->outIdx;
        e1->side = EdgeSide::Right;
        e2->side = EdgeSide::Left;
        e = e2;
        if(e->prevInAEL == e1)
            prevE = e1->prevInAEL;
        else prevE = e->prevInAEL;
    }

    if(prevE && prevE->outIdx >= 0 && math::LT(prevE->top[1], pt[1]) && math::LT(e->top[1], pt[1])) {
        auto xPrev = TopX(*prevE, pt[1]);
        auto xE = TopX(*e, pt[1]);
        if(math::EQ(xPrev, xE) && (e->windDelta != 0) && (prevE->windDelta != 0) &&
            SlopesEqual(Point<num_type>(xPrev, pt[1]), prevE->top, Point<num_type>(xE, pt[1]), e->top, m_useFullRange)) {
            auto * outPt = AddOutPt(prevE, pt);
            AddJoin(result, outPt, e->top);
        }
    }
    return result;
}

template <typename num_type>
inline OutRec<num_type> * Clipper<num_type>::GetOutRec(size_t idx)
{
    auto * outrec = m_polyOuts[idx];
    while(outrec != m_polyOuts[outrec->idx])
        outrec = m_polyOuts[outrec->idx];
    return outrec;
}

template <typename num_type>
inline void Clipper<num_type>::AppendPolygon(TEdge<num_type> * e1, TEdge<num_type> * e2)
{
    //get the start and ends of both output polygons ...
    auto * outRec1 = m_polyOuts[e1->outIdx];
    auto * outRec2 = m_polyOuts[e2->outIdx];

    OutRec<num_type> * holeStateRec = nullptr;
    if(OutRec1RightOfOutRec2(outRec1, outRec2))
        holeStateRec = outRec2;
    else if(OutRec1RightOfOutRec2(outRec2, outRec1))
        holeStateRec = outRec1;
    else
        holeStateRec = GetLowermostRec(outRec1, outRec2);

    //get the start and ends of both output polygons and
    //join e2 poly onto e1 poly and delete pointers to e2 ...

    auto * p1Lft = outRec1->pts;
    auto * p1Rt = p1Lft->prev;
    auto * p2Lft = outRec2->pts;
    auto * p2Rt = p2Lft->prev;

    //join e2 poly onto e1 poly and delete pointers to e2 ...
    if(e1->side == EdgeSide::Left) {
        if(e2->side == EdgeSide::Left) {
            //z y x a b c
            ReversePolyPtLinks(p2Lft);
            p2Lft->next = p1Lft;
            p1Lft->prev = p2Lft;
            p1Rt->next = p2Rt;
            p2Rt->prev = p1Rt;
            outRec1->pts = p2Rt;
        }
        else {
            //x y z a b c
            p2Rt->next = p1Lft;
            p1Lft->prev = p2Rt;
            p2Lft->prev = p1Rt;
            p1Rt->next = p2Lft;
            outRec1->pts = p2Lft;
        }
    }
    else {
        if(e2->side == EdgeSide::Right) {
            //a b c z y x
            ReversePolyPtLinks(p2Lft);
            p1Rt->next = p2Rt;
            p2Rt->prev = p1Rt;
            p2Lft->next = p1Lft;
            p1Lft->prev = p2Lft;
        }
        else {
            //a b c x y z
            p1Rt->next = p2Lft;
            p2Lft->prev = p1Rt;
            p1Lft->prev = p2Rt;
            p2Rt->next = p1Lft;
        }
    }

    outRec1->bottomPt = nullptr;
    if(holeStateRec == outRec2) {
        if(outRec2->firstLeft != outRec1)
            outRec1->firstLeft = outRec2->firstLeft;
        outRec1->isHole = outRec2->isHole;
    }
    outRec2->pts = nullptr;
    outRec2->bottomPt = nullptr;
    outRec2->firstLeft = outRec1;

    int oKIdx = e1->outIdx;
    int obsoleteIdx = e2->outIdx;

    e1->outIdx = constant::unassigned; //nb: safe because we only get here via AddLocalMaxPoly
    e2->outIdx = constant::unassigned;

    auto * e = m_activeEdges;
    while(e) {
        if(e->outIdx == obsoleteIdx) {
            e->outIdx = oKIdx;
            e->side = e1->side;
            break;
        }
        e = e->nextInAEL;
    }

    outRec2->idx = outRec1->idx;
}

template <typename num_type>
inline void Clipper<num_type>::IntersectEdges(TEdge<num_type> * e1, TEdge<num_type> * e2, Point<num_type> & pt)
{
    bool e1Contributing = e1->outIdx >= 0;
    bool e2Contributing = e2->outIdx >= 0;

    //if either edge is on an OPEN path ...
    if (e1->windDelta == 0 || e2->windDelta == 0) {
        //ignore subject-subject open path intersections UNLESS they
        //are both open paths, AND they are both 'contributing maximas' ...
        if(e1->windDelta == 0 && e2->windDelta == 0) return;

        //if intersecting a subj line with a subj poly ...
        else if(e1->polyTyp == e2->polyTyp && e1->windDelta != e2->windDelta && m_clipType == ClipType::Union) {
            if (e1->windDelta == 0) {
                if (e2Contributing) {
                    AddOutPt(e1, pt);
                    if (e1Contributing)
                        e1->outIdx = constant::unassigned;
                }
            }
            else {
                if (e1Contributing) {
                    AddOutPt(e2, pt);
                    if(e2Contributing)
                        e2->outIdx = constant::unassigned;
                }
            }
        }
        else if(e1->polyTyp != e2->polyTyp) {
            //toggle subj open path OutIdx on/off when Abs(clip.WndCnt) == 1 ...
            if ((e1->windDelta == 0) && std::abs(e2->windCnt) == 1 &&
                    (m_clipType != ClipType::Union || e2->windCnt2 == 0)) {
                AddOutPt(e1, pt);
                if(e1Contributing)
                    e1->outIdx = constant::unassigned;
            }
            else if(e2->windDelta == 0 && std::abs(e1->windCnt) == 1 &&
                    (m_clipType != ClipType::Union || e1->windCnt2 == 0)) {
                AddOutPt(e2, pt);
                if (e2Contributing)
                    e2->outIdx = constant::unassigned;
            }
        }
        return;
    }

    //update winding counts...
    //assumes that e1 will be to the Right of e2 ABOVE the intersection
    if(e1->polyTyp == e2->polyTyp) {
        if(isEvenOddFillType( *e1)) {
            int oldE1WindCnt = e1->windCnt;
            e1->windCnt = e2->windCnt;
            e2->windCnt = oldE1WindCnt;
        }
        else {
            if(e1->windCnt + e2->windDelta == 0)
                e1->windCnt = -e1->windCnt;
            else e1->windCnt += e2->windDelta;
            if(e2->windCnt - e1->windDelta == 0)
                e2->windCnt = -e2->windCnt;
            else e2->windCnt -= e1->windDelta;
        }
    }
    else {
        if(!isEvenOddFillType(*e2))
            e1->windCnt2 += e2->windDelta;
        else e1->windCnt2 = e1->windCnt2 == 0 ? 1 : 0;
        if(!isEvenOddFillType(*e1))
            e2->windCnt2 -= e1->windDelta;
        else e2->windCnt2 = e2->windCnt2 == 0 ? 1 : 0;
    }

    PolyFillType e1FillType1, e2FillType1, e1FillType2, e2FillType2;
    if (e1->polyTyp == PolyType::Subject) {
        e1FillType1 = m_subjFillType;
        e1FillType2 = m_clipFillType;
    }
    else {
        e1FillType1 = m_clipFillType;
        e1FillType2 = m_subjFillType;
    }
    if(e2->polyTyp == PolyType::Subject) {
        e2FillType1 = m_subjFillType;
        e2FillType2 = m_clipFillType;
    }
    else {
        e2FillType1 = m_clipFillType;
        e2FillType2 = m_subjFillType;
    }

    int e1Wc, e2Wc;
    switch(e1FillType1) {
        case PolyFillType::Positive : 
            e1Wc =  e1->windCnt; break;
        case PolyFillType::Negative :
            e1Wc = -e1->windCnt; break;
        default :
            e1Wc = std::abs(e1->windCnt);
    }
    switch(e2FillType1) {
        case PolyFillType::Positive :
            e2Wc =  e2->windCnt; break;
        case PolyFillType::Negative :
            e2Wc = -e2->windCnt; break;
        default :
            e2Wc = std::abs(e2->windCnt);
    }

    if(e1Contributing && e2Contributing) {
        if((e1Wc != 0 && e1Wc != 1) ||
            (e2Wc != 0 && e2Wc != 1) ||
            (e1->polyTyp != e2->polyTyp && m_clipType != ClipType::Xor)) {
            AddLocalMaxPoly(e1, e2, pt);
        }
        else {
            AddOutPt(e1, pt);
            AddOutPt(e2, pt);
            SwapSides(*e1 ,*e2);
            SwapPolyIndexes(*e1 ,*e2);
        }
    }
    else if(e1Contributing) {
        if(e2Wc == 0 || e2Wc == 1) {
            AddOutPt(e1, pt);
            SwapSides(*e1, *e2);
            SwapPolyIndexes(*e1, *e2);
        }
    }
    else if(e2Contributing) {
        if(e1Wc == 0 || e1Wc == 1) {
            AddOutPt(e2, pt);
            SwapSides(*e1, *e2);
            SwapPolyIndexes(*e1, *e2);
        }
    }
    else if((e1Wc == 0 || e1Wc == 1) && (e2Wc == 0 || e2Wc == 1)) {
        //neither edge is currently contributing ...
        int e1Wc2, e2Wc2;
        switch(e1FillType2) {
            case PolyFillType::Positive :
                e1Wc2 =  e1->windCnt2; break;
            case PolyFillType::Negative :
                e1Wc2 = -e1->windCnt2; break;
            default :
                e1Wc2 = std::abs(e1->windCnt2);
        }
        switch(e2FillType2) {
            case PolyFillType::Positive :
                e2Wc2 =  e2->windCnt2; break;
            case PolyFillType::Negative :
                e2Wc2 = -e2->windCnt2; break;
            default :
                e2Wc2 = std::abs(e2->windCnt2);
        }
        if(e1->polyTyp != e2->polyTyp) {
            AddLocalMinPoly(e1, e2, pt);
        }
        else if(e1Wc == 1 && e2Wc == 1) {
            switch(m_clipType) {
                case ClipType::Intersection :
                    if(e1Wc2 > 0 && e2Wc2 > 0)
                        AddLocalMinPoly(e1, e2, pt);
                    break;
                case ClipType::Union :
                    if(e1Wc2 <= 0 && e2Wc2 <= 0)
                        AddLocalMinPoly(e1, e2, pt);
                    break;
                case ClipType::Difference :
                    if(((e1->polyTyp == PolyType::Clip) && (e1Wc2 > 0) && (e2Wc2 > 0)) ||
                        ((e1->polyTyp == PolyType::Subject) && (e1Wc2 <= 0) && (e2Wc2 <= 0)))
                            AddLocalMinPoly(e1, e2, pt);
                    break;
                case ClipType::Xor :
                    AddLocalMinPoly(e1, e2, pt);
                    break;
            }
        }
        else SwapSides(*e1, *e2);
    }    
}

template <typename num_type>
inline OutPt<num_type> * Clipper<num_type>::AddOutPt(TEdge<num_type> * e, const Point<num_type> & pt)
{
    if(e->outIdx < 0) {
        auto * outRec = ClipperBase<num_type>::CreateOutRec();
        outRec->isOpen = e->windDelta == 0;
        auto * newOp = new OutPt<num_type>;
        outRec->pts = newOp;
        newOp->idx = outRec->idx;
        newOp->pt = pt;
        newOp->next = newOp;
        newOp->prev = newOp;
        if(!outRec->isOpen)
            SetHoleState(e, outRec);
        e->outIdx = outRec->idx;
        return newOp;
    }
    else {
        auto * outRec = m_polyOuts[e->outIdx];
        //OutRec.pts is the 'Left-most' point & OutRec.pts.prev is the 'Right-most'
        auto * op = outRec->pts;

        bool toFront = e->side == EdgeSide::Left;
        if(toFront && pt == op->pt) return op;
        else if(!toFront && pt == op->prev->pt) return op->prev;

        auto * newOp = new OutPt<num_type>;
        newOp->idx = outRec->idx;
        newOp->pt = pt;
        newOp->next = op;
        newOp->prev = op->prev;
        newOp->prev->next = newOp;
        op->prev = newOp;
        if(toFront) outRec->pts = newOp;
        return newOp;
    }
}

template <typename num_type>
inline OutPt<num_type> * Clipper<num_type>::GetLastOutPt(TEdge<num_type> * e)
{
    auto * outRec = m_polyOuts[e->outIdx];
    if(e->side == EdgeSide::Left) return outRec->pts;
    else return outRec->pts->prev;
}

template <typename num_type>
inline bool Clipper<num_type>::ProcessIntersections(const num_type topY)
{
    if(!m_activeEdges ) return true;
    try {
        BuildIntersectList(topY);
        auto ilSize = m_intersectList.size();
        if(ilSize == 0) return true;
        if(ilSize == 1 || FixupIntersectionOrder())
            ProcessIntersectList();
        else return false;
    }
    catch(...) {
        m_sortedEdges = nullptr;
        DisposeIntersectNodes();
        throw Exception("ProcessIntersections error");
    }
    m_sortedEdges = nullptr;
    return true;   
}

template <typename num_type>
inline void Clipper<num_type>::SetHoleState(TEdge<num_type> * e, OutRec<num_type> * outrec)
{
    auto * e2 = e->prevInAEL;
    TEdge<num_type> * eTmp = nullptr;
    while(e2) {
        if (e2->outIdx >= 0 && e2->windDelta != 0) {
            if(!eTmp) eTmp = e2;
            else if (eTmp->outIdx == e2->outIdx) eTmp = 0;
        }
        e2 = e2->prevInAEL;
    }
    if(!eTmp) {
        outrec->firstLeft = nullptr;
        outrec->isHole = false;
    }
    else {
        outrec->firstLeft = m_polyOuts[eTmp->outIdx];
        outrec->isHole = !(outrec->firstLeft->isHole);
    } 
}

template <typename num_type>
inline void Clipper<num_type>::DisposeIntersectNodes()
{
    for(auto * node : m_intersectList)
        delete node;
    m_intersectList.clear();
}

template <typename num_type>
inline void Clipper<num_type>::FixupOutPolygon(OutRec<num_type> & outrec)
{
    //FixupOutPolygon() - removes duplicate points and simplifies consecutive
    //parallel edges by removing the middle vertex.
    OutPt<num_type> * lastOK = nullptr;
    outrec.bottomPt = nullptr;
    auto * pp = outrec.pts;
    bool preserveCol = m_preserveCollinear || m_strictSimple;

    for(;;) {
        if(pp->prev == pp || pp->prev == pp->next) {
            DisposeOutPts(pp);
            outrec.pts = nullptr;
            return;
        }

        //test for duplicate points and collinear edges ...
        if(pp->pt == pp->next->pt || pp->pt == pp->prev->pt ||
            (SlopesEqual(pp->prev->pt, pp->pt, pp->next->pt, m_useFullRange) &&
            (!preserveCol || !Pt2IsBetweenPt1AndPt3(pp->prev->pt, pp->pt, pp->next->pt)))) {
            lastOK = nullptr;
            auto * tmp = pp;
            pp->prev->next = pp->next;
            pp->next->prev = pp->prev;
            pp = pp->prev;
            delete tmp;
        }
        else if(pp == lastOK) break;
        else {
            if(!lastOK) lastOK = pp;
            pp = pp->next;
        }
    }
    outrec.pts = pp;
}
    
template <typename num_type>
inline void Clipper<num_type>::FixupOutPolyline(OutRec<num_type> & outrec)
{
    auto * pp = outrec.pts;
    auto * lastPP = pp->prev;
    while(pp != lastPP) {
        pp = pp->next;
        if(pp->pt == pp->prev->pt) {
            if(pp == lastPP) lastPP = pp->prev;
            auto * tmpPP = pp->prev;
            tmpPP->next = pp->next;
            pp->next->prev = tmpPP;
            delete pp;
            pp = tmpPP;
        }
    }

    if(pp == pp->prev) {
        DisposeOutPts(pp);
        outrec.pts = nullptr;
        return;
    }
}

template <typename num_type>
inline bool Clipper<num_type>::FixupIntersectionOrder()
{
    //pre-condition: intersections are sorted Bottom-most first.
    //Now it's crucial that intersections are made only between adjacent edges,
    //so to ensure this the order of intersections may need adjusting ...
    CopyAELToSEL();
    using Node = IntersectNode<num_type>;
    auto cmp = [](Node * n1, Node * n2){ return math::LT(n2->point[1], n1->point[1]); };
    std::sort(m_intersectList.begin(), m_intersectList.end(), cmp);
    size_t size = m_intersectList.size();
    for (size_t i = 0; i < size; ++i) {
        if(!EdgesAdjacent(*m_intersectList[i])) {
            size_t j = i + 1;
            while (j < size && !EdgesAdjacent(*m_intersectList[j])) j++;
            if (j == size)  return false;
            std::swap(m_intersectList[i], m_intersectList[j]);
        }
        SwapPositionsInSEL(m_intersectList[i]->edge1, m_intersectList[i]->edge2);
    }
    return true;
}

template <typename num_type>
inline void Clipper<num_type>::FixHoleLinkage(OutRec<num_type> & outrec)
{
    //skip OutRecs that (a) contain outermost polygons or
    //(b) already have the correct owner/child linkage ...
    if(!outrec.firstLeft || (outrec.isHole != outrec.firstLeft->isHole && outrec.firstLeft->pts)) return;

    auto * orfl = outrec.firstLeft;
    while(orfl && ((orfl->isHole == outrec.isHole) || !orfl->pts))
        orfl = orfl->firstLeft;
    outrec.firstLeft = orfl;
}

template <typename num_type>
inline void Clipper<num_type>::AddJoin(OutPt<num_type> * op1, OutPt<num_type> * op2, Point<num_type> offPt)
{
    auto * j = new Join<num_type>;
    j->outPt1 = op1;
    j->outPt2 = op2;
    j->offPt = std::move(offPt);
    m_joins.push_back(j);
}

template <typename num_type>
inline void Clipper<num_type>::AddGhostJoin(OutPt<num_type> * op, Point<num_type> offPt)
{
    auto * j = new Join<num_type>;
    j->outPt1 = op;
    j->outPt2 = 0;
    j->offPt = std::move(offPt);
    m_ghostJoins.push_back(j);
}

template <typename num_type>
inline void Clipper<num_type>::ClearJoins()
{
    for(auto * join : m_joins)
        delete join;
    m_joins.clear();
}

template <typename num_type>
inline void Clipper<num_type>::ClearGhostJoins()
{
    for(auto * ghostJoin : m_ghostJoins)
        delete ghostJoin;
    m_ghostJoins.clear();
}

template <typename num_type>
inline bool Clipper<num_type>::JoinPoints(Join<num_type> * j, OutRec<num_type> * outRec1, OutRec<num_type> * outRec2)
{
    OutPt<num_type> * op1 = j->outPt1, * op1b = nullptr;
    OutPt<num_type> * op2 = j->outPt2, * op2b = nullptr;

    //There are 3 kinds of joins for output polygons ...
    //1. Horizontal joins where Join.OutPt1 & Join.OutPt2 are vertices anywhere
    //along (horizontal) collinear edges (& Join.OffPt is on the same horizontal).
    //2. Non-horizontal joins where Join.OutPt1 & Join.OutPt2 are at the same
    //location at the Bottom of the overlapping segment (& Join.OffPt is above).
    //3. StrictSimple joins where edges touch but are not collinear and where
    //Join.OutPt1, Join.OutPt2 & Join.OffPt all share the same point.
    bool isHorizontal = math::EQ(j->outPt1->pt[1], j->offPt[1]);

    if(isHorizontal && (j->offPt == j->outPt1->pt) && (j->offPt == j->outPt2->pt)) {
        //Strictly Simple join ...
        if(outRec1 != outRec2) return false;
        op1b = j->outPt1->next;
        while(op1b != op1 && (op1b->pt == j->offPt)) op1b = op1b->next;
        bool reverse1 = math::GT(op1b->pt[1], j->offPt[1]);
        op2b = j->outPt2->next;
        while(op2b != op2 && (op2b->pt == j->offPt)) op2b = op2b->next;
        bool reverse2 = math::GT(op2b->pt[1], j->offPt[1]);
        if(reverse1 == reverse2) return false;
        if(reverse1) {
            op1b = DupOutPt(op1, false);
            op2b = DupOutPt(op2, true );
            op1->prev = op2;
            op2->next = op1;
            op1b->next = op2b;
            op2b->prev = op1b;
            j->outPt1 = op1;
            j->outPt2 = op1b;
            return true;
        }
        else {
            op1b = DupOutPt(op1, true );
            op2b = DupOutPt(op2, false);
            op1->next = op2;
            op2->prev = op1;
            op1b->prev = op2b;
            op2b->next = op1b;
            j->outPt1 = op1;
            j->outPt2 = op1b;
            return true;
        }
    }
    else if(isHorizontal) {
        //treat horizontal joins differently to non-horizontal joins since with
        //them we're not yet sure where the overlapping is. OutPt1.Pt & OutPt2.Pt
        //may be anywhere along the horizontal edge.
        op1b = op1;
        while(op1->prev->pt[1] == op1->pt[1] && op1->prev != op1b && op1->prev != op2) op1 = op1->prev;
        while(op1b->next->pt[1] == op1b->pt[1] && op1b->next != op1 && op1b->next != op2) op1b = op1b->next;
        if(op1b->next == op1 || op1b->next == op2) return false; //a flat 'polygon'

        op2b = op2;
        while(op2->prev->pt[1] == op2->pt[1] && op2->prev != op2b && op2->prev != op1b) op2 = op2->prev;
        while(op2b->next->pt[1] == op2b->pt[1] && op2b->next != op2 && op2b->next != op1) op2b = op2b->next;
        if(op2b->next == op2 || op2b->next == op1) return false; //a flat 'polygon'

        num_type left, right;
        //Op1 --> Op1b & Op2 --> Op2b are the extremites of the horizontal edges
        if(!GetOverlap(op1->pt[0], op1b->pt[0], op2->pt[0], op2b->pt[0], left, right)) return false;

        //DiscardLeftSide: when overlapping edges are joined, a spike will created
        //which needs to be cleaned up. However, we don't want Op1 or Op2 caught up
        //on the discard Side as either may still be needed for other joins ...
        Point<num_type> pt;
        bool discardLeftSide;
        if(math::GE(op1->pt[0], left) && math::LE(op1->pt[0], right)) {
            pt = op1->pt; discardLeftSide = (op1->pt[0] > op1b->pt[0]);
        }
        else if(math::GE(op2->pt[0], left) && math::LE(op2->pt[0], right)) {
            pt = op2->pt; discardLeftSide = (op2->pt[0] > op2b->pt[0]);
        }
        else if(math::GE(op1b->pt[0], left) && math::LE(op1b->pt[0], right)) {
            pt = op1b->pt; discardLeftSide = op1b->pt[0] > op1->pt[0];
        }
        else {
            pt = op2b->pt; discardLeftSide = math::GT(op2b->pt[0], op2->pt[0]);
        }
        j->outPt1 = op1; j->outPt2 = op2;
        return JoinHorz(op1, op1b, op2, op2b, pt, discardLeftSide);
    }
    else {
        //nb: For non-horizontal joins ...
        //1. Jr.OutPt1.Pt[1] == Jr.OutPt2.Pt[1]
        //2. Jr.OutPt1.Pt > Jr.OffPt[1]

        //make sure the polygons are correctly oriented ...
        op1b = op1->next;
        while(op1b->pt == op1->pt && op1b != op1) op1b = op1b->next;
        bool reverse1 = math::GT(op1b->pt[1], op1->pt[1]) || !SlopesEqual(op1->pt, op1b->pt, j->offPt, m_useFullRange);
        if(reverse1) {
            op1b = op1->prev;
            while((op1b->pt == op1->pt) && (op1b != op1)) op1b = op1b->prev;
            if(math::GT(op1b->pt[1], op1->pt[1]) || !SlopesEqual(op1->pt, op1b->pt, j->offPt, m_useFullRange)) return false;
        }
        op2b = op2->next;
        while(op2b->pt == op2->pt && op2b != op2) op2b = op2b->next;
        bool reverse2 = math::GT(op2b->pt[1], op2->pt[1]) || !SlopesEqual(op2->pt, op2b->pt, j->offPt, m_useFullRange);
        if(reverse2) {
            op2b = op2->prev;
            while(op2b->pt == op2->pt && op2b != op2) op2b = op2b->prev;
            if(math::GT(op2b->pt[1], op2->pt[1]) || !SlopesEqual(op2->pt, op2b->pt, j->offPt, m_useFullRange)) return false;
        }

        if(op1b == op1 || op2b == op2 || op1b == op2b || (outRec1 == outRec2 && reverse1 == reverse2)) return false;

        if(reverse1) {
            op1b = DupOutPt(op1, false);
            op2b = DupOutPt(op2, true );
            op1->prev = op2;
            op2->next = op1;
            op1b->next = op2b;
            op2b->prev = op1b;
            j->outPt1 = op1;
            j->outPt2 = op1b;
            return true;
        }
        else {
            op1b = DupOutPt(op1, true );
            op2b = DupOutPt(op2, false);
            op1->next = op2;
            op2->prev = op1;
            op1b->prev = op2b;
            op2b->next = op1b;
            j->outPt1 = op1;
            j->outPt2 = op1b;
            return true;
        }
    } 
}

template <typename num_type>
inline void Clipper<num_type>::JoinCommonEdges()
{
    for (auto * join : m_joins) {
        auto * outRec1 = GetOutRec(join->outPt1->idx);
        auto * outRec2 = GetOutRec(join->outPt2->idx);

        if(!outRec1->pts || !outRec2->pts) continue;
        if(outRec1->isOpen || outRec2->isOpen) continue;

        //get the polygon fragment with the correct hole state (FirstLeft)
        //before calling JoinPoints() ...
        OutRec<num_type> * holeStateRec = nullptr;
        if(outRec1 == outRec2) holeStateRec = outRec1;
        else if(OutRec1RightOfOutRec2(outRec1, outRec2)) holeStateRec = outRec2;
        else if(OutRec1RightOfOutRec2(outRec2, outRec1)) holeStateRec = outRec1;
        else holeStateRec = GetLowermostRec(outRec1, outRec2);

        if(!JoinPoints(join, outRec1, outRec2)) continue;

        if(outRec1 == outRec2) {
            //instead of joining two polygons, we've just created a new one by
            //splitting one polygon into two.
            outRec1->pts = join->outPt1;
            outRec1->bottomPt = nullptr;
            outRec2 = ClipperBase<num_type>::CreateOutRec();
            outRec2->pts = join->outPt2;

            //update all OutRec2.Pts Idx's ...
            UpdateOutPtIdxs(*outRec2);

            if(Poly2ContainsPoly1(outRec2->pts, outRec1->pts)) {
                //outRec1 contains outRec2 ...
                outRec2->isHole = !outRec1->isHole;
                outRec2->firstLeft = outRec1;

                if(m_usingPolyTree) FixupFirstLefts2(outRec2, outRec1);

                if((outRec2->isHole ^ m_reverseOutput) == math::isPositive(Area(*outRec2)))
                    ReversePolyPtLinks(outRec2->pts);

            }
            else if(Poly2ContainsPoly1(outRec1->pts, outRec2->pts)) {
                //outRec2 contains outRec1 ...
                outRec2->isHole = outRec1->isHole;
                outRec1->isHole = !outRec2->isHole;
                outRec2->firstLeft = outRec1->firstLeft;
                outRec1->firstLeft = outRec2;

                if(m_usingPolyTree) FixupFirstLefts2(outRec1, outRec2);

                if((outRec1->isHole ^ m_reverseOutput) == math::isPositive(Area(*outRec1)))
                    ReversePolyPtLinks(outRec1->pts);
            }
            else {
                //the 2 polygons are completely separate ...
                outRec2->isHole = outRec1->isHole;
                outRec2->firstLeft = outRec1->firstLeft;

                //fixup FirstLeft pointers that may need reassigning to OutRec2
                if(m_usingPolyTree) FixupFirstLefts1(outRec1, outRec2);
            }

        }
        else {
            //joined 2 polygons together ...
            outRec2->pts = nullptr;
            outRec2->bottomPt = nullptr;
            outRec2->idx = outRec1->idx;

            outRec1->isHole = holeStateRec->isHole;
            if(holeStateRec == outRec2)
                outRec1->firstLeft = outRec2->firstLeft;
            outRec2->firstLeft = outRec1;

            if(m_usingPolyTree) FixupFirstLefts3(outRec2, outRec1);
        }
    }
}

template <typename num_type>
inline void Clipper<num_type>::DoSimplePolygons()
{
    size_t i = 0;
    while(i < m_polyOuts.size()){
        auto * outrec = m_polyOuts[i++];
        auto * op = outrec->pts;
        if(!op || outrec->isOpen) continue;
        do {//for each Pt in Polygon until duplicate found do ...
            auto * op2 = op->next;
            while(op2 != outrec->pts) {
                if(op->pt == op2->pt && op2->next != op && op2->prev != op) {
                    //split the polygon into two ...
                    auto * op3 = op->prev;
                    auto * op4 = op2->prev;
                    op->prev = op4;
                    op4->next = op;
                    op2->prev = op3;
                    op3->next = op2;

                    outrec->pts = op;
                    auto * outrec2 = ClipperBase<num_type>::CreateOutRec();
                    outrec2->pts = op2;
                    UpdateOutPtIdxs(*outrec2);
                    if(Poly2ContainsPoly1(outrec2->pts, outrec->pts)) {
                        //OutRec2 is contained by OutRec1 ...
                        outrec2->isHole = !outrec->isHole;
                        outrec2->firstLeft = outrec;
                        if(m_usingPolyTree) FixupFirstLefts2(outrec2, outrec);
                    }
                    else if(Poly2ContainsPoly1(outrec->pts, outrec2->pts)) {
                        //OutRec1 is contained by OutRec2 ...
                        outrec2->isHole = outrec->isHole;
                        outrec->isHole = !outrec2->isHole;
                        outrec2->firstLeft = outrec->firstLeft;
                        outrec->firstLeft = outrec2;
                        if(m_usingPolyTree) FixupFirstLefts2(outrec, outrec2);
                    }
                    else {
                        //the 2 polygons are separate ...
                        outrec2->isHole = outrec->isHole;
                        outrec2->firstLeft = outrec->firstLeft;
                        if(m_usingPolyTree) FixupFirstLefts1(outrec, outrec2);
                    }
                    op2 = op;//ie get ready for the Next iteration
                }
                op2 = op2->next;
            }
            op = op->next;
        } while (op != outrec->pts);
    }
}

template <typename num_type>
inline void Clipper<num_type>::FixupFirstLefts1(OutRec<num_type> * oldOutRec, OutRec<num_type> * newOutRec)
{
    //tests if NewOutRec contains the polygon before reassigning FirstLeft
    for(auto * outRec : m_polyOuts) {
        auto * firstLeft = ParseFirstLeft(outRec->firstLeft);
        if(outRec->pts && firstLeft == oldOutRec) {
        if(Poly2ContainsPoly1(outRec->pts, newOutRec->pts))
            outRec->firstLeft = newOutRec;
        }
    }
}

template <typename num_type>
inline void Clipper<num_type>::FixupFirstLefts2(OutRec<num_type> * innerOutRec, OutRec<num_type> * outerOutRec)
{
    //A polygon has split into two such that one is now the inner of the other.
    //It's possible that these polygons now wrap around other polygons, so check
    //every polygon that's also contained by OuterOutRec's FirstLeft container
    //(including 0) to see if they've become inner to the new inner polygon ...
    auto * orfl = outerOutRec->firstLeft;
    for (auto * outRec : m_polyOuts) {
        if(!outRec->pts || outRec == outerOutRec || outRec == innerOutRec) continue;
        auto * firstLeft = ParseFirstLeft(outRec->firstLeft);
        if(firstLeft != orfl && firstLeft != innerOutRec && firstLeft != outerOutRec) continue;
        if(Poly2ContainsPoly1(outRec->pts, innerOutRec->pts))
            outRec->firstLeft = innerOutRec;
        else if(Poly2ContainsPoly1(outRec->pts, outerOutRec->pts))
            outRec->firstLeft = outerOutRec;
        else if(outRec->firstLeft == innerOutRec || outRec->firstLeft == outerOutRec)
            outRec->firstLeft = orfl;
    }
}

template <typename num_type>
inline void Clipper<num_type>::FixupFirstLefts3(OutRec<num_type> * oldOutRec, OutRec<num_type> * newOutRec)
{
    //reassigns FirstLeft WITHOUT testing if NewOutRec contains the polygon
    for(auto * outRec : m_polyOuts) {
        auto * firstLeft = ParseFirstLeft(outRec->firstLeft);
        if(outRec->pts && firstLeft == oldOutRec)
            outRec->firstLeft = newOutRec;
    }
}

template <typename num_type>
inline void Clipper<num_type>::InsertEdgeIntoAEL(TEdge<num_type> * edge, TEdge<num_type> * startEdge)
{
    if(!m_activeEdges) {
        edge->prevInAEL = nullptr;
        edge->nextInAEL = nullptr;
        m_activeEdges = edge;
    }
    else if(!startEdge && E2InsertsBeforeE1(*m_activeEdges, *edge)) {
        edge->prevInAEL = nullptr;
        edge->nextInAEL = m_activeEdges;
        m_activeEdges->prevInAEL = edge;
        m_activeEdges = edge;
    }
    else {
        if(!startEdge) startEdge = m_activeEdges;
        while(startEdge->nextInAEL && !E2InsertsBeforeE1(*startEdge->nextInAEL , *edge)) {
            startEdge = startEdge->nextInAEL;
        }
        edge->nextInAEL = startEdge->nextInAEL;
        if(startEdge->nextInAEL)
            startEdge->nextInAEL->prevInAEL = edge;
        edge->prevInAEL = startEdge;
        startEdge->nextInAEL = edge;
    }
}

template <typename num_type>
inline void Clipper<num_type>::AddEdgeToSEL(TEdge<num_type> * edge)
{
    //SEL pointers in PEdge are reused to build a list of horizontal edges.
    //However, we don't need to worry about order with horizontal edge processing.
    if(!m_sortedEdges) {
        m_sortedEdges = edge;
        edge->prevInSEL = nullptr;
        edge->nextInSEL = nullptr;
    }
    else {
        edge->nextInSEL = m_sortedEdges;
        edge->prevInSEL = nullptr;
        m_sortedEdges->prevInSEL = edge;
        m_sortedEdges = edge;
    }
}

template <typename num_type>
inline bool Clipper<num_type>::PopEdgeFromSEL(TEdge<num_type> * & edge)
{
    if(!m_sortedEdges) return false;
    edge = m_sortedEdges;
    DeleteFromSEL(m_sortedEdges);
    return true;
}

template <typename num_type>
inline void Clipper<num_type>::CopyAELToSEL()
{
    auto * e = m_activeEdges;
    m_sortedEdges = e;
    while(e) {
        e->prevInSEL = e->prevInAEL;
        e->nextInSEL = e->nextInAEL;
        e = e->nextInAEL;
    }
}

template <typename num_type>
inline void Clipper<num_type>::DeleteFromSEL(TEdge<num_type> * e)
{
    auto * selPrev = e->prevInSEL;
    auto * selNext = e->nextInSEL;
    if(!selPrev && !selNext && e != m_sortedEdges) return; //already deleted
    if(selPrev)
        selPrev->nextInSEL = selNext;
    else m_sortedEdges = selNext;
    if(selNext)
        selNext->prevInSEL = selPrev;
    e->nextInSEL = nullptr;
    e->prevInSEL = nullptr;
}

}//namespace clipper
}//namespace geometry
}//namespace generic