/**
 * @file PolygonMerge.hpp
 * @author bwu
 * @brief Utility for polygon merge
 * @version 0.1
 * @date 2022-02-25
 */
#ifndef GENERIC_GEOMETRY_POLYGONMERGE_HPP
#define GENERIC_GEOMETRY_POLYGONMERGE_HPP
#define GENERIC_GEOMETRY_POLYGONMERGE_USE_RTREE
#include "generic/thread/TaskFlow.hpp"
#include "BooleanOperation.hpp"
#include "Connectivity.hpp"
#include "HashFunction.hpp"
#include "Utility.hpp"

#ifdef  GENERIC_GEOMETRY_POLYGONMERGE_USE_RTREE
#include "generic/tree/RectTree.hpp"
#else
#include "generic/tree/QuadTreeUtilityMT.hpp"
#endif//GENERIC_GEOMETRY_POLYGONMERGE_USE_RTREE

namespace generic {
namespace geometry{

template <typename property_type, typename num_type>
class PolygonWithProp
{
public:
    using float_t = common::float_type<num_type>;
    property_type property;
    Polygon2D<num_type> solid;
    std::list<Polygon2D<num_type> > holes;

    bool hasHole() const
    {
        return holes.size() > 0;
    }

    Box2D<num_type> BBox() const
    {
        Box2D<num_type> bbox;
        bbox |= Extent(solid);
        for(const auto & hole : holes)
            bbox |= Extent(hole);
        return bbox;
    }

    void Normalize()
    {
        if(!solid.isCCW())
            solid.Reverse();
        for(auto & hole : holes)
            if(hole.isCCW()) hole.Reverse();
    }

    float_t CoveredArea() const
    {
        return boost::polygon::area(solid);
    }

    void RemoveTinyHoles(float_t area)
    {
        for(auto iter = holes.begin(); iter != holes.end();) {
            if(math::LT<float_t>(boost::polygon::area(*iter), area)) {
                iter = holes.erase(iter);
            }
            else iter++;
        }
    }
};

template <typename property_type, typename num_type>
struct PolygonWithPropExt
{
    Box2D<num_type> operator() (const PolygonWithProp<property_type, num_type> & p) const
    {
        return p.BBox();
    }
};

#ifdef GENERIC_GEOMETRY_POLYGONMERGE_USE_RTREE
template <typename property_type, typename num_type>
using PolygonMergeTaskTree = tree::RectTree<num_type, PolygonWithProp<property_type, num_type>, PolygonWithPropExt<property_type, num_type> >;

template <typename property_type, typename num_type>
using PolygonMergeTaskNode = tree::RectNode<num_type, PolygonWithProp<property_type, num_type> >;

template <typename property_type, typename num_type>
using PolygonMergeSubTaskNodes = typename tree::RectNode<num_type, PolygonWithProp<property_type, num_type> >::RectChildren;

#else
template <typename property_type, typename num_type>
using PolygonMergeTaskTree = tree::QuadTree<num_type, PolygonWithProp<property_type, num_type>, PolygonWithPropExt<property_type, num_type> >;

template <typename property_type, typename num_type>
using PolygonMergeTaskNode = tree::QuadTree<num_type, PolygonWithProp<property_type, num_type>, PolygonWithPropExt<property_type, num_type> >;

template <typename property_type, typename num_type>
using PolygonMergeSubTaskNodes = typename tree::QuadTree<num_type, PolygonWithProp<property_type, num_type>, PolygonWithPropExt<property_type, num_type> >::QuadChildren;

#endif//GENERIC_GEOMETRY_POLYGONMERGE_USE_RTREE

template <typename property_type, typename num_type>
using PropDiffPolygon = std::pair<std::set<property_type>, std::list<Polygon2D<num_type> > >;
template <typename property_type, typename num_type>
using PropDiffPolygons = std::list<PropDiffPolygon<property_type, num_type> >;

namespace detail {

template <typename num_type, typename = void>
class Point2DIndexMap {};

template <typename num_type>
class Point2DIndexMap<num_type, typename std::enable_if<std::is_integral<num_type>::value>::type>
{
public:
    using Point = Point2D<num_type>;

    void Clear() { m_ptMap.clear(); }
    void Insert(const Point & p, size_t index) { m_ptMap.insert(std::make_pair(p, index)); }
    bool Find(const Point & p, size_t & index) const
    {
        for(int i : {0, -1, 1}) {
            for(int j : {0, -1, 1}) {
                auto iter = m_ptMap.find(Point(p[0] + i, p[1] + j));
                if(iter != m_ptMap.end()) {
                    index = iter->second;
                    return true;
                }
            }
        }
        return false;
    }
private:
    std::unordered_map<Point, size_t, PointHash<num_type> > m_ptMap;
};

template <typename num_type>
class Point2DIndexMap<num_type, typename std::enable_if<std::is_floating_point<num_type>::value>::type>
{
public:
    using Point = Point2D<num_type>;
    using IndexPoint = std::pair<Point, size_t>;

    void Clear() { m_rtree.clear(); }
    void Insert(const Point & p, size_t index) { m_rtree.insert(std::make_pair(p, index)); }
    bool Find(const Point & p, size_t & index) const
    {
        std::vector<IndexPoint> query;
        m_rtree.query(boost::geometry::index::nearest(p, 1), std::back_inserter(query));
        if(query.empty()) return false;
        if(query.front().first == p) {
            index = query.front().second;
            return true;
        }
        return false;
    }
private:
    boost::geometry::index::rtree<IndexPoint, boost::geometry::index::linear<16, 4> > m_rtree;
};

}//namespace detail

template <typename Merger> class PolygonMergeRunner;
template <typename property_type, typename num_type>
class PolygonMerger
{
    friend class PolygonMergeRunner<PolygonMerger<property_type, num_type> >;
public:
    using float_t = common::float_type<num_type>;
    using PolygonData = PolygonWithProp<property_type, num_type>;
    using PropDiffAreas = PropDiffPolygons<property_type, num_type>;
    using MergeTaskTree = PolygonMergeTaskTree<property_type, num_type>;
    using MergeTaskNode = PolygonMergeTaskNode<property_type, num_type>;
    using MergeSubTaskNodes = PolygonMergeSubTaskNodes<property_type, num_type>;
    using PropertyMap = std::unordered_map<property_type, property_type>;

    struct MergeSettings
    {
        bool cleanPolyonPoints = false;
        bool checkPropertyDiff = false;
        bool ignoreTinySolid = false;
        bool ignoreTinyHoles = false;
        float_t tinySolidArea = 0;
        float_t tinyHolesArea = 0;
        float_t cleanPointDist = 0;
        size_t mergeThreashold = 1024;
    };
    ~PolygonMerger() { Clear(); }

    void SetMergeSettings(const MergeSettings & settings);

    PolygonData * AddObject(property_type property, Box2D<num_type> box);
    PolygonData * AddObject(property_type property, Polygon2D<num_type> polygon);
    PolygonData * AddObject(property_type property, PolygonWithHoles2D<num_type> pwh);

    void Merge();//single thread

    void GetAllPolygons(std::list<PolygonData * > & polygons);

    void GetAllPolygons(std::list<const PolygonData * > & polygons);

    const Box2D<num_type> & GetBBox() const { return m_bbox; }

    const PropDiffAreas & GetPropDiffAreas() const { return m_propDiffAreas; }

    void Clear();
private:    
    void PreProcess();
    MergeTaskTree & GetMergeTaskTree() { return m_mergeTaskTree; }
    void MergeRegion(MergeTaskNode * node);
    void PostProcess();

private:
    void BuildTaskTree();
    void CleanPolygons();
    void FilterOutTinyArea();
    void FilterOutTinyArea(std::list<PolygonData * > & polygons);
    void FilterOutTinyHoles(std::list<PolygonData * > & polygons);
    void GetOverlappedSubTaskNodes(MergeTaskNode * node, std::vector<std::list<MergeTaskNode * > > & nodeGroups);
    void ExtractAndMergeComponentsPolygons(std::list<PolygonData * > & polygons);
    void ExtractComponentsByPolygonPropertye(const std::vector<PolygonData * > & vecPolys, std::vector<std::vector<size_t> > & cc);
    void ExtractComponentsByPhysicConnection(const std::vector<PolygonData * > & vecPolys, std::vector<std::vector<size_t> > & cc);
    void MergeComponentsPolygons(std::list<PolygonData * > & polygons, const std::vector<PolygonData * > & vecPolys, const std::vector<std::vector<size_t> > & cc);
    void MergePolygons(std::list<PolygonData * > & polygons);

    PolygonData * AddPolygonData(PolygonData * pd);
    PolygonData * makePolygonData(Polygon2D<num_type> & in, property_type prop);
    static void Simplify(Polygon2D<num_type> & polygon, std::list<Polygon2D<num_type> > & holes);

private:
    Box2D<num_type> m_bbox;
    std::list<PolygonData * > m_datas;

    PropertyMap   m_propertyMap;
    MergeTaskTree m_mergeTaskTree;
    MergeSettings m_mergeSettings;
    PropDiffAreas m_propDiffAreas;
};

template <typename property_type, typename num_type>
class PolygonMergeUtils
{
public:
    using float_t = typename PolygonMerger<property_type, num_type>::float_t;
    using PolygonData = typename PolygonMerger<property_type, num_type>::PolygonData;
    static void CleanPolygons(std::list<PolygonData * > & polygons, const float_t dist);
    static void CleanPolygon(PolygonData * polygon, const float_t dist);
    static void CleanPolygon(Polygon2D<num_type> & polygon, const float_t dist);
};

template <typename Merger>
class PolygonMergeRunner
{
    using TaskNode = thread::taskflow::TaskNode;
    using TaskFlow = thread::taskflow::TaskFlow;
    using MergeTaskNode = typename Merger::MergeTaskNode;
    using MergeTaskTree = typename Merger::MergeTaskTree;

public:
    PolygonMergeRunner(Merger & merger, size_t threads)
     :m_merger(merger), m_threads(threads) {}
    ~PolygonMergeRunner() = default;

    void Run();
private:
    void ScheduleTasks(MergeTaskTree & tree);
    void ScheduleSubTasks(MergeTaskNode * parent, TaskNode * successor);

private:
    std::unique_ptr<TaskFlow> m_taskFlow = nullptr;
    Merger & m_merger;
    size_t m_threads;
};

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::SetMergeSettings(const MergeSettings & settings)
{
    m_mergeSettings = settings;
}

template <typename property_type, typename num_type>
inline typename PolygonMerger<property_type, num_type>::PolygonData *
PolygonMerger<property_type, num_type>::AddObject(property_type property, Box2D<num_type> box)
{
    auto data = new PolygonData;
    data->property = property;
    data->solid = toPolygon(box);
    return AddPolygonData(data);
}

template <typename property_type, typename num_type>
inline typename PolygonMerger<property_type, num_type>::PolygonData *
PolygonMerger<property_type, num_type>::AddObject(property_type property, Polygon2D<num_type> polygon)
{
    auto data = new PolygonData;
    data->property = property;
    data->solid = std::move(polygon);
    return AddPolygonData(data);
}

template <typename property_type, typename num_type>
inline typename PolygonMerger<property_type, num_type>::PolygonData *
PolygonMerger<property_type, num_type>::AddObject(property_type property, PolygonWithHoles2D<num_type> pwh)
{
    auto data = new PolygonData;
    data->property = property;
    data->solid = std::move(pwh.outline);
    for(auto & hole : pwh.holes)
        data->holes.emplace_back(std::move(hole));

    return AddPolygonData(data); 
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::Merge()//single thread
{
    PreProcess();
    MergeRegion(&m_mergeTaskTree);
    PostProcess();
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::GetAllPolygons(std::list<PolygonData * > & polygons)
{
    polygons.clear();
    m_mergeTaskTree.GetAllObjects(polygons);
    if(polygons.empty())
        polygons.insert(polygons.end(), m_datas.begin(), m_datas.end());
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::GetAllPolygons(std::list<const PolygonData * > & polygons)
{
    polygons.clear();
    m_mergeTaskTree.GetAllObjects(polygons);
    if(polygons.empty())
        polygons.insert(polygons.end(), m_datas.begin(), m_datas.end());
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::PreProcess()
{
    if(m_mergeSettings.cleanPolyonPoints &&
        math::isPositive(m_mergeSettings.cleanPointDist))
        CleanPolygons();
    
    BuildTaskTree();
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::MergeRegion(MergeTaskNode * node)
{
    if(node->hasChild()) {
        const auto & children = node->GetChildren();
        for(auto * child : children)
            MergeRegion(child);
    }

    bool merged = false;
    std::list<PolygonData * > allObjs;
    std::list<PolygonData * > mergedObjs;
#ifdef GENERIC_GEOMETRY_POLYGONMERGE_USE_RTREE
    std::vector<std::list<MergeTaskNode * > > subNodeGroups;
    GetOverlappedSubTaskNodes(node, subNodeGroups);
    if(!subNodeGroups.empty()) {
        for(const auto & subNodeGroup : subNodeGroups) {
            std::list<PolygonData * > objs;
            for(auto * subNode : subNodeGroup) {
                const auto & subObjs = subNode->GetObjs();
                objs.insert(objs.end(), subObjs.begin(), subObjs.end());
                subNode->Clear();
            }
            ExtractAndMergeComponentsPolygons(objs);
            mergedObjs.insert(mergedObjs.end(), objs.begin(), objs.end());
        }
        merged = true;
    }
#endif//GENERIC_GEOMETRY_POLYGONMERGE_USE_RTREE
    node->GetAllObjects(allObjs);

    if(mergedObjs.size())
        allObjs.insert(allObjs.end(), mergedObjs.begin(), mergedObjs.end());
    
    if(false == node->GetObjs().empty()) {
        ExtractAndMergeComponentsPolygons(allObjs);
        merged = true;
    }

    if(merged) FilterOutTinyHoles(allObjs);
    node->Build(std::move(allObjs));
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::PostProcess()
{
    if(m_mergeSettings.cleanPointDist && math::isPositive(m_mergeSettings.cleanPointDist)) CleanPolygons();
    if(m_mergeSettings.ignoreTinySolid && math::isPositive(m_mergeSettings.tinySolidArea)) FilterOutTinyArea();
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::BuildTaskTree()
{
    m_mergeTaskTree.SetBBox(m_bbox);
    m_mergeTaskTree.Build(m_datas, m_mergeSettings.mergeThreashold);
    m_datas.clear();
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::CleanPolygons()
{
    std::list<PolygonData * > polygons;
    GetAllPolygons(polygons);
    PolygonMergeUtils<property_type, num_type>::CleanPolygons(polygons, m_mergeSettings.cleanPointDist);
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::FilterOutTinyArea()
{
    std::list<PolygonData * > polygons;
    GetAllPolygons(polygons);
    FilterOutTinyArea(polygons);
    m_mergeTaskTree.Build(polygons, 0);
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::FilterOutTinyArea(std::list<PolygonData * > & polygons)
{
    if(m_mergeSettings.ignoreTinySolid && math::isPositive(m_mergeSettings.tinySolidArea)) {
        for(auto iter = polygons.begin(); iter != polygons.end();) {
            if(math::LT((*iter)->CoveredArea(), m_mergeSettings.tinySolidArea)) {
                delete (*iter);
                iter = polygons.erase(iter);
            }
            else iter++;
        }
    }
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::FilterOutTinyHoles(std::list<PolygonData * > & polygons)
{
    if(m_mergeSettings.ignoreTinyHoles && math::isPositive(m_mergeSettings.tinyHolesArea)) {
        for(auto * polygon : polygons)
            polygon->RemoveTinyHoles(m_mergeSettings.tinyHolesArea);
    }
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::GetOverlappedSubTaskNodes(MergeTaskNode * node, std::vector<std::list<MergeTaskNode * > > & nodeGroups)
{
    using namespace topology;
    auto subNodes = node->GetChildren();
    std::vector<MergeTaskNode * > vecNodes;
    vecNodes.reserve(subNodes.size());
    for(auto * node : subNodes) vecNodes.push_back(node);
    
    std::vector<std::set<int> > connection;
    auto boxGetter = [](MergeTaskNode * node) { return node->GetBBox(); };
    ConnectivityExtraction(vecNodes, boxGetter, connection);

    auto graph = makeSparseIndexGraph(connection);

    std::vector<std::list<size_t> > cc;
    ConnectedComponents(*graph, cc);

    for(const auto & component : cc) {
        if(component.size() > 1) {
            nodeGroups.push_back(std::list<MergeTaskNode * >{});
            for(auto i : component) {
                nodeGroups.back().push_back(vecNodes[i]);
            }
        }
    }
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::ExtractAndMergeComponentsPolygons(std::list<PolygonData * > & polygons)
{
    if(polygons.size() <= 1) return;

    std::vector<PolygonData *> vecPolys;
    vecPolys.assign(polygons.begin(), polygons.end());

    std::vector<std::vector<size_t> > cc;
    if(m_mergeSettings.checkPropertyDiff) {
        ExtractComponentsByPhysicConnection(vecPolys, cc);
    }
    else {
        ExtractComponentsByPolygonPropertye(vecPolys, cc);
    }

    MergeComponentsPolygons(polygons, vecPolys, cc);
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::ExtractComponentsByPolygonPropertye(const std::vector<PolygonData * > & vecPolys, std::vector<std::vector<size_t> > & cc)
{
    cc.clear();
    std::unordered_map<property_type, size_t> ccMap;
    for(size_t i = 0; i < vecPolys.size(); ++i) {
        const auto & prop = vecPolys[i]->property;
        auto iter = ccMap.find(prop);
        if(iter == ccMap.end()){
            iter = ccMap.insert(std::make_pair(prop, ccMap.size())).first;
            cc.push_back(std::vector<size_t>{});
        }
        cc[iter->second].push_back(i);
    }
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::ExtractComponentsByPhysicConnection(const std::vector<PolygonData * > & vecPolys, std::vector<std::vector<size_t> > & cc)
{
    auto dup = [](const Polygon2D<num_type> & from, Polygon2D<num_type> & to) {
        std::copy(from.GetPoints().begin(), from.GetPoints().end(), std::back_inserter(to.GetPoints()));
    };

    auto geomGetter = [&dup](const PolygonData * pd) {
        PolygonWithHoles2D<num_type> pwh;
        dup(pd->solid, pwh.outline);
        for(const auto & hole : pd->holes) {
            pwh.holes.push_back(Polygon2D<num_type>{});
            dup(hole, pwh.holes.back());
        }
        return pwh;
    };

    using namespace topology;
    std::vector<std::set<int> > connection;
    ConnectivityExtraction(vecPolys, geomGetter, connection);

    auto g = topology::makeSparseIndexGraph(connection);

    ConnectedComponents(*g, cc);
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::MergeComponentsPolygons(std::list<PolygonData * > & polygons, const std::vector<PolygonData * > & vecPolys, const std::vector<std::vector<size_t> > & cc)
{
    polygons.clear();
    std::list<PolygonData * > tmp;
    for(const auto & c : cc) {
        if(c.size() == 1) polygons.push_back(vecPolys[c.front()]);
        else {
            tmp.clear();
            for(auto i : c)
                tmp.push_back(vecPolys[i]);
            
            MergePolygons(tmp);
            polygons.insert(polygons.end(), tmp.begin(), tmp.end());
        }
    }
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::MergePolygons(std::list<PolygonData * > & polygons)
{
    if(polygons.size() <= 1) return;
    boost::polygon::property_merge<num_type, property_type> merger;
    std::map<std::set<property_type>, boolean::PolygonSet2D<num_type> > results;
    for(auto * pd : polygons) {
        auto property = m_propertyMap.count(pd->property) ?
                        m_propertyMap.at(pd->property) : pd->property;
        merger.insert(pd->solid, property, false);
        for(const auto & hole : pd->holes) {
            merger.insert(hole, property, true);
        }
        delete pd;
    }

    polygons.clear();
    merger.merge(results);

    for(auto & result : results) {
        std::list<Polygon2D<num_type> > outs;

        if(!result.second.empty())
            result.second.get(outs);

        auto & properties = result.first;
        GENERIC_ASSERT(!properties.empty())
        if(properties.size() > 1) {
            if(m_mergeSettings.checkPropertyDiff) {
                m_propDiffAreas.emplace_back(std::make_pair(std::move(properties), std::move(outs)));
                outs.clear();
            }
            else {
                auto iter = properties.begin();
                auto prop = *iter; iter++;
                for(; iter != properties.end(); ++iter) {
                    m_propertyMap.insert(std::make_pair(*iter, prop));
                }
            }
        }
        auto prop = *(properties.begin());
        for(auto & out : outs) {
            polygons.push_back(makePolygonData(out, prop));
        }
    }
}

template <typename property_type, typename num_type>
inline typename PolygonMerger<property_type, num_type>::PolygonData * 
PolygonMerger<property_type, num_type>::AddPolygonData(PolygonData * pd)
{
    pd->Normalize();
    m_bbox |= pd->BBox();
    m_datas.push_back(pd);
    return pd;
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::Simplify(Polygon2D<num_type> & polygon, std::list<Polygon2D<num_type> > & holes)
{
    holes.clear();
    if(polygon.Front() == polygon.Back()) polygon.PopBack();
    if(polygon.Size() <= 3) return;

    size_t size = polygon.Size();
    struct PtNode { size_t prev; size_t next; };
    std::vector<PtNode> nodeList(size);
    for(size_t i = 0; i < size; ++i) {
        nodeList[i].prev = (i + size - 1) % size;
        nodeList[i].next = (i + 1) % size;
    }

    size_t res;
    std::list<Polygon2D<num_type> > tmp;
    detail::Point2DIndexMap<num_type> ptMap;
    for(size_t i = 0; i < polygon.Size(); ++i) {
        auto found = ptMap.Find(polygon[i], res);
        if(found) {
            size_t prev = res;
            size_t curr = i;
            size_t next = nodeList[curr].next;

            nodeList[curr].next = prev;
            Polygon2D<num_type> hole;
            size_t start = prev;
            size_t index = start;
            while(start != nodeList[index].next) {
                hole << polygon[index];
                index = nodeList[index].next;
            }
            Simplify(hole, tmp);
            if(hole.Size() >= 3)
                holes.emplace_back(std::move(hole));
            prev = nodeList[prev].prev;
            nodeList[prev].next = curr;
            nodeList[curr].prev = prev;
            nodeList[curr].next = next;
        }
        ptMap.Insert(polygon[i], i);
    }

    size_t start = size - 1;
    size_t index = start;
    Polygon2D<num_type> simplified;
    while(start != nodeList[index].next) {
        simplified << polygon[index];
        index = nodeList[index].next;
    }
    simplified << polygon[index];
    std::swap(polygon, simplified);
}

template <typename property_type, typename num_type>
inline typename PolygonMerger<property_type, num_type>::PolygonData * 
PolygonMerger<property_type, num_type>::makePolygonData(Polygon2D<num_type> & in, property_type prop)
{
    auto pd = new PolygonData;
    pd->property = prop;
    pd->solid = std::move(in);
    Simplify(pd->solid, pd->holes);
    return pd;
}

template <typename property_type, typename num_type>
inline void PolygonMerger<property_type, num_type>::Clear()
{
    std::list<PolygonData * > polygons;
    GetAllPolygons(polygons);
    PolygonData * pd = nullptr;
    auto iter = polygons.begin();
    for(; iter != polygons.end(); ++iter) {
        pd = *iter;
        if(pd) delete pd;
    }

    m_datas.clear();
    m_propertyMap.clear();
    m_propDiffAreas.clear();
    m_mergeTaskTree.Clear();
    m_bbox.SetInvalid();
}
template <typename property_type, typename num_type>
inline void PolygonMergeUtils<property_type, num_type>::CleanPolygons(std::list<PolygonData * > & polygons, float_t dist)
{
    for(auto & polygon : polygons)
        CleanPolygon(polygon, dist);
}

template <typename property_type, typename num_type>
inline void PolygonMergeUtils<property_type, num_type>::CleanPolygon(PolygonData * polygon, const float_t dist)
{
    CleanPolygon(polygon->solid, dist);
    for(auto & hole : polygon->holes)
        CleanPolygon(hole, dist);
}

template <typename property_type, typename num_type>
inline void PolygonMergeUtils<property_type, num_type>::CleanPolygon(Polygon2D<num_type> & polygon, const float_t dist)
{
    Polygon2D<num_type> in, out;
    out = polygon;
    size_t size = out.Size();
    while(size != in.Size()){
        in = std::move(out);
        boost::geometry::simplify(in, out, dist);
        size = out.Size();
    }
    if(DistanceSq(out.Front(), out.Back()) > dist * dist) out.PopBack();
    if(out.Size() >= 3) polygon = std::move(out);
}

template <typename Merger>
inline void PolygonMergeRunner<Merger>::Run()
{
    if(m_threads > 1) {
        m_merger.PreProcess();

        auto & topNode = m_merger.GetMergeTaskTree();
        ScheduleTasks(topNode);

        thread::taskflow::Executor executor(m_threads);
        executor.Run(*m_taskFlow);

        m_merger.PostProcess();
    }
    else {
        m_merger.Merge();
    }
}

template <typename Merger>
inline void PolygonMergeRunner<Merger>::ScheduleTasks(MergeTaskTree & tree)
{
    m_taskFlow.reset(new TaskFlow);
    auto * task = m_taskFlow->Emplace(std::bind(&Merger::MergeRegion, &m_merger, &tree));
    ScheduleSubTasks(&tree, task);
}

template <typename Merger>
inline void PolygonMergeRunner<Merger>::ScheduleSubTasks(MergeTaskNode * parent, TaskNode * successor)
{
    if(!(parent->hasChild())) return;
    const auto & children = parent->GetChildren();
    for(auto * child : children){
        auto * task = m_taskFlow->Emplace(std::bind(&Merger::MergeRegion, &m_merger, child));
        task->Precede(successor);
        ScheduleSubTasks(child, task);
    }
}

}//namespace geometry
}//namespace generic
#endif//GENERIC_GEOMETRY_POLYGONMERGE_HPP
