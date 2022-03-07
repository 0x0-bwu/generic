/**
 * @file Clipper.hpp
 * @author bwu
 * @brief Clipper library
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_GEOMETRY_CLIPPER_CLIPPER_HPP
#define GENERIC_GEOMETRY_CLIPPER_CLIPPER_HPP
#include <boost/multiprecision/cpp_int.hpp>
#include "generic/common/Exception.hpp"
#include "Polygon.hpp"

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

enum class ClipType { Intersection, Union, Difference, Xor };
enum class EndType  { ClosedPolygon, ClosedLine, OpenButt, OpenSquare, OpenRound };
enum class JoinType { Square, Round, Miter };
enum class PolyType { Subject, Clip };
enum class PolyFillType { EvenOdd, NonZero, Positive, Negative };
enum class EdgeSide { Left = 1, Right = 2 };

template <typename num_type>
using Point = Point2D<num_type>;

template <typename num_type>
using Path = Polyline2D<num_type>;

template <typename num_type>
using Paths = std::vector<Path<num_type> >;

template <typename num_type>
class PolyNode
{
public:
    using Children = std::vector<PolyNode<num_type> * >;
    PolyNode() = default;
    virtual ~PolyNode() {};
    Path<num_type> contour;
    Children children;
    PolyNode<num_type> * parent = nullptr;
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
    child.index = index;
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
  bool isHole;
  bool isOpen;
  OutRec<num_type> * firstLeft = nullptr;  //see comments in clipper.pas
  PolyNode<num_type> * polyNd = nullptr;
  OutPt<num_type> * pts = nullptr;
  OutPt<num_type> * bottomPt = nullptr;
};

template <typename num_type>
struct TEdge
{
    using float_t = float_type<num_type>;
    Point<num_type> bot;
    Point<num_type> curr; //current (updated for every new scanbeam)
    Point<num_type> Top;
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
    TEdge<num_type> * nrevInAEL = nullptr;
    TEdge<num_type> * nextInSEL = nullptr;
    TEdge<num_type> * prevInSEL = nullptr;
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
    num_type dy = (e.top[1] - e.Bot[1]);
    if(dy == 0) e.dx = constant::horizontal;
    else e.dx = static_cast<float_t>(e.top[0] - e.bot[0]) / dy;
}

template <typename num_type>
inline void InitEdge2(TEdge<num_type> & e, PolyType polyTyp)
{
    if(e.curr[1] >= e.next->curr[1]) {
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
    std::swap(e.top[0], e.bot[1]);
}

template <typename num_type>
inline bool isHorizontal(TEdge<num_type> & e)
{
    using float_t = float_type<num_type>;
    return math::EQ<float_t>(e.dx, constant::horizontal);
}

template <typename num_type>
inline TEdge<num_type> * FindNextLocMin(TEdge<num_type> * e)
{
    for(;;) {
        while (e->bot != e->prev->bot || e->curr == e->top) e = e->next;
        if (!isHorizontal(*e) && !isHorizontal(*e->prev)) break;
        while (isHorizontal(*e->prev)) e = e->prev;
        auto * e2 = e;
        while (isHorizontal(*e)) e = e->next;
        if (e->top[1] == e->prev->bot[1]) continue; //ie just an intermediate horz.
        if (e2->prev->bot[0] < e->bot[0]) e = e2;
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
struct LocalMinimum {
    num_type y;
    TEdge<num_type> * boundL = nullptr;
    TEdge<num_type> * boundR = nullptr;
};

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
inline bool SlopesEqual(const Point<num_type> p1, const Point<num_type> p2, const Point<num_type> p3, bool useFullInt64Range)
{
    if constexpr (std::is_same<num_type, int64_t>::value) {
        if(useFullInt64Range)
            return Int128Mul(p1[1] - p2[1], p2[0] - p3[0]) == Int128Mul(p1[0] - p2[0], p2[1] - p3[1]);
        else return (p1[1] - p2[1]) * (p2[0] - p3[0]) == (p1[0] - p2[0]) * (p2[1] - p3[1]);
    }
    else return (p1[1] - p2[1]) * (p2[0] - p3[0]) == (p1[0] - p2[0]) * (p2[1] - p3[1]);
}

template <typename num_type>
inline bool Pt2IsBetweenPt1AndPt3(const Point<num_type> & p1, const Point<num_type> & p2, const Point<num_type> & p3)
{
  if((p1 == p3) || (p1 == p2) || (p3 == p2)) return false;
  else if (p1[0] != p3[0]) return (p2[0] > p1[0]) == (p2[0] < p3[0]);
  else return (p2[1] > p1[1]) == (p2[1] < p3[1]);
}

template <typename num_type>
using EdgeList = std::vector<TEdge<num_type> * >;

template <typename num_type>
using PolyOutList = std::vector<OutRec<num_type> * >;

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
    // virtual void Reset();
    TEdge<num_type> * ProcessBound(TEdge<num_type> * e, bool IsClockwise);
    // void InsertScanbeam(const cInt Y);
    // bool PopScanbeam(cInt &Y);
    // bool LocalMinimaPending();
    // bool PopLocalMinima(cInt Y, const LocalMinimum *&locMin);
    // OutRec* CreateOutRec();
    // void DisposeAllOutRecs();
    // void DisposeOutRec(PolyOutList::size_type index);
    // void SwapPositionsInAEL(TEdge *edge1, TEdge *edge2);
    // void DeleteFromAEL(TEdge *e);
    // void UpdateEdgeIntoAEL(TEdge *&e);

protected:
    MinimaListIter m_currentLM;
    MinimaList m_minimaList;
    bool m_useFullRange;
    EdgeList<num_type> m_edges;
    bool m_preserveCollinear;
    bool m_hasOpenPaths;
    PolyOutList<num_type> m_polyOuts;
    TEdge<num_type> * m_activeEdges;
    ScanbeamList     m_scanbeam;
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
inline bool ClipperBase<num_type>::AddPath(const Path<num_type> & path, PolyType polyTyp, bool closed)
{
    if(!closed && polyTyp == PolyType::Clip)
        throw Exception("AddPath: Open paths must be subject.");

    int highI = static_cast<int>(path.size()) - 1;
    if (closed) while (highI > 0 && (path[highI] == path[0])) --highI;
    while (highI > 0 && (path[highI] == path[highI -1])) --highI;
    if ((closed && highI < 2) || (!closed && highI < 1)) return false;

    //create a new edge array ...
    auto * edges = new TEdge[highI +1];

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
            if(e == e->Next) break;
            if(e == eStart) eStart = e->next;
            e = RemoveEdge(e);
            eLoopStop = e;
            continue;
        }
        if(e->prev == e->next) break; //only two vertices
        else if(closed && SlopesEqual(e->rrev->curr, e->curr, e->next->curr, m_useFullRange) &&
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
        locMin[1] = e->bot[1];
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
        locMin[1] = e->bot[1];
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

        e = ProcessBound(locMin.LeftBound, leftBoundIsForward);
        if (e->OutIdx == constant::skip) e = ProcessBound(e, leftBoundIsForward);

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
            while (e != result && isHorizontal(*e)) e = e->Prev;
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
            locMin[1] = e->bot[1];
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
            e = e->Next;
        }
        if(isHorizontal(*e) && e != eStart && e->bot[0] != e->prev->top[0])
            ReverseHorizontal(*e);
        result = result->next; //move to the edge just beyond current bound
    }
    else {
        while(result->top[1] == result->prev->Bot[1] && result->prev->outIdx != constant::skip)
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

}//namespace clipper
}//namespace geometry
}//namespace generic

#endif //GENERIC_GEOMETRY_CLIPPER_CLIPPER_HPP