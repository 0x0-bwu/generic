#ifndef GENERIC_TREE_QUADTREE_HPP
#define GENERIC_TREE_QUADTREE_HPP
#include "generic/common/Traits.hpp"
#include "generic/geometry/Box.hpp"
#include "generic/geometry/Utility.hpp"
#include <array>
#include <list>
#include <map>
namespace generic{
namespace tree {
using generic::geometry::Box2D;
using generic::geometry::Point2D;
using generic::common::num_integer_tag;
using generic::common::num_floating_tag;
using generic::common::num_traits_float_or_int;
template <typename num_type, typename object, typename extent>
class QuadTree
{
public:
    enum class Direction { Sourth = 0, East = 1, North = 2, West = 3};
    enum class Orientation {SW = 0, SE = 1, NE = 2, NW = 3};
    using QuadNode = QuadTree<num_type, object, extent>;
    using QuadChildren = std::array<QuadNode * , 4>;
    QuadTree():m_level(0), m_parent(nullptr){ m_children = {nullptr}; }
    QuadTree(const Box2D<num_type> & bbox, QuadNode * parent = nullptr): m_bbox(bbox), m_parent(parent)
    {
        m_children = {nullptr};
        m_level = parent ? parent->GetLevel() + 1 : 0;
    }
    ~QuadTree(){ Clear(); }
    
    void Build(size_t max_objs_to_build_sub);
    void Build(std::list<object* > && objs);
    void Build(std::list<object* > objs, size_t max_objs_to_build_sub);
    bool Insert(object * obj, extent & ext);
    bool CreateSubNodes(size_t max_objs_to_build_sub);
    bool CreateSubNodes();

    bool FillObjects(std::list<object * > & objs);

    int GetLevel() const { return m_level; }
    void SetBBox(const Box2D<num_type > & bbox) { m_bbox = bbox; }
    const Box2D<num_type > & GetBBox() const { return m_bbox; }

    std::list<object * > & GetObjs() { return m_objs; }
    const std::list<object * > & GetObjs() const { return m_objs; }

    bool isRoot() const { return nullptr == m_parent; }
    QuadNode * GetParent() { return m_parent; }
    bool hasChild() const { return m_children.front() != nullptr; }
    QuadNode * GetChild(Orientation o) { return hasChild() ? m_children[ static_cast<size_t>(o) ] : nullptr; }
    QuadChildren & GetChildren() { return m_children; }
    const QuadChildren & GetChildren() const { return m_children; }
    void GetAllObjects(std::list<object * > & allObjs) const;
    void GetAllObjects(std::list<const object * > & allObjs) const;

    void Split();
    void Balance();

    static void GetAllLeafNodes(QuadNode * node, std::list<QuadNode * > & allNodes);
    static void GetAllNodesByLevel(QuadNode * node, std::map<int, std::list<QuadNode * > > & allNodes);
    static QuadNode * Neighbor(QuadNode * node, QuadNode * root, Direction direction);
    
    template <typename Condition>
    static void CreateSubNodesIf(QuadNode * node, Condition && condition);

    template <typename T>
    bool CreateSubBoxCondition() const;
    bool CreateSubBoxCondition(const num_floating_tag) const;
    bool CreateSubBoxCondition(const num_integer_tag) const;

    void ClearSubNodes();
    void Clear();

private:
    int m_level;
    Box2D<num_type > m_bbox;
    std::list<object * > m_objs;
    QuadNode * m_parent;
    QuadChildren m_children;
};

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent >::ClearSubNodes()
{
    for(size_t i = 0; i < 4; ++i){
        if(m_children[ i ]) {
            delete m_children[ i ];
            m_children[ i ] = nullptr;
        }
    }
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent >::Clear()
{
    ClearSubNodes();
    m_objs.clear();
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::Build(size_t max_objs_to_build_sub)
{
    if(CreateSubNodes(max_objs_to_build_sub)){
        for(auto * child : m_children)
            child->Build(max_objs_to_build_sub);
    }
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::Build(std::list<object* > && objs)
{
    Clear();
    extent ext;
    m_objs = objs;
    for(auto * obj : m_objs)
        m_bbox |= ext(*obj);
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::Build(std::list<object * > objs, size_t max_objs_to_build_sub)
{
    Build(std::move(objs));
    Build(max_objs_to_build_sub);
}

template <typename num_type, typename object, typename extent>
inline bool QuadTree<num_type, object, extent>::Insert(object * obj, extent & ext)
{
    if(m_bbox >= ext(*obj)){
        m_objs.push_back(obj);
        return true;
    }
    return false;
}

template <typename num_type, typename object, typename extent>
inline bool QuadTree<num_type, object, extent>::CreateSubNodes(size_t max_objs_to_build_sub)
{
    if(max_objs_to_build_sub == 0) return false;
    if(max_objs_to_build_sub >= m_objs.size()) return false;
    if(!CreateSubBoxCondition<num_type>()) return false;

    return CreateSubNodes();
}

template <typename num_type, typename object, typename extent>
inline bool QuadTree<num_type, object, extent>::CreateSubNodes()
{
    Split();
    bool res = false;
    for(auto * child : m_children)
        res |= child->FillObjects(m_objs);
    if(false == res){
        ClearSubNodes();
        return false;
    }
    return true;
}

template <typename num_type, typename object, typename extent>
inline bool QuadTree<num_type, object, extent>::FillObjects(std::list<object * > & objs)
{
    extent ext;
    bool res = false;
    auto iter = objs.begin();
    while(iter != objs.end()){
        if(this->Insert(*iter, ext)){
            iter = objs.erase(iter);
            res = true;
        }
        else iter++;
    }
    return res;
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::GetAllObjects(std::list<object * > & allObjs) const
{
    for(auto * obj : m_objs)
        allObjs.push_back(obj);

    if(hasChild()){
        for(auto * child : m_children)
            child->GetAllObjects(allObjs);
    }
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::GetAllObjects(std::list<const object * > & allObjs) const
{
    for(auto * obj : m_objs)
        allObjs.push_back(obj);

    if(hasChild()){
        for(auto * child : m_children)
            child->GetAllObjects(allObjs);
    }
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::Split()
{
    ClearSubNodes();
    Box2D<num_type> bbox[4];
    Point2D<num_type> ct = m_bbox.Center().template Cast<num_type>();
    Point2D<num_type> ll = Point2D<num_type>(m_bbox[0][0], m_bbox[0][1]);
    Point2D<num_type> lm = Point2D<num_type>(ct[0], m_bbox[0][1]);
    Point2D<num_type> mr = Point2D<num_type>(m_bbox[1][0], ct[1]);
    Point2D<num_type> ur = Point2D<num_type>(m_bbox[1][0], m_bbox[1][1]);
    Point2D<num_type> um = Point2D<num_type>(ct[0], m_bbox[1][1]);
    Point2D<num_type> ml = Point2D<num_type>(m_bbox[0][0], ct[1]);
    bbox[0] = Box2D<num_type>(ll, ct);
    bbox[1] = Box2D<num_type>(lm, mr);
    bbox[2] = Box2D<num_type>(ct, ur);
    bbox[3] = Box2D<num_type>(ml, um);

    for(int i = 0; i < 4; ++i){
        m_children[i] = new QuadTree<num_type, object, extent>(bbox[i], this);
    }
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::Balance()
{
    std::list<QuadNode * > nodeList;
    GetAllLeafNodes(this, nodeList);

    typename std::list<QuadNode * >::iterator iter = nodeList.begin();
    while(nodeList.size()){
        QuadNode * node = *iter;

        for(size_t i = 0; i < 4; ++i){
            Direction d = static_cast<Direction>(i);
            QuadNode * neighbor = Neighbor(node, this, d);
            if(neighbor && neighbor->hasChild()){
                Orientation o1, o2;
                if(Direction::Sourth == d){ o1 = Orientation::NW; o2 = Orientation::NE; }
                else if(Direction::East == d){ o1 = Orientation::NW; o2 = Orientation::SW; }
                else if(Direction::North == d){ o1 = Orientation::SW; o2 = Orientation::SE; }
                else { o1 = Orientation::NE; o2 = Orientation::SE; }
                if(neighbor->GetChild(o1)->hasChild() || neighbor->GetChild(o2)->hasChild()){
                    node->Split();
                    const QuadChildren & children = node->GetChildren();
                    for(size_t i = 0; i < children.size(); ++i){
                        nodeList.push_back(children[ i ]);
                    }

                    if(node->GetLevel() > neighbor->GetLevel()){
                        neighbor->Split();
                        const QuadChildren & children = neighbor->GetChildren();
                        for(size_t i = 0; i < children.size(); ++i){
                            nodeList.push_back(children[ i ]);
                        }
                    }
                }
            }
        }
        iter = nodeList.erase(iter);
    }
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::GetAllLeafNodes(QuadNode * node, std::list<QuadNode * > & allLeaves)
{
    if(!node->hasChild()) {
        allLeaves.push_back(node);
        return;
    }
    const QuadChildren & children = node->GetChildren();
    for(auto * child : children) GetAllLeafNodes(child, allLeaves);
}

template <typename num_type, typename object, typename extent>
inline void QuadTree<num_type, object, extent>::GetAllNodesByLevel(QuadNode * node, std::map<int, std::list<QuadNode * > >& allNodes)
{
    int topLevel = node->GetLevel();
    allNodes[topLevel].push_back(node);
    if(!node->hasChild()) return;
    const QuadChildren & children = node->GetChildren();
    for(auto * child : children) GetAllNodesByLevel(child, allNodes);
}

template <typename num_type, typename object, typename extent>
inline QuadTree<num_type, object, extent > * QuadTree<num_type, object, extent>::Neighbor(QuadNode * node, QuadNode * root, Direction direction)
{
    if(node == root) return nullptr;
    QuadNode * parent = node->GetParent();
    if(Direction::Sourth == direction){
        if(node == parent->GetChild(Orientation::NW)) return parent->GetChild(Orientation::SW);
        if(node == parent->GetChild(Orientation::NE)) return parent->GetChild(Orientation::SE);

        QuadNode * pare_neighbor = Neighbor(parent, root, Direction::Sourth);
        if(nullptr == pare_neighbor || !pare_neighbor->hasChild()) return pare_neighbor;
        else{
            if(node == parent->GetChild(Orientation::SW)) return pare_neighbor->GetChild(Orientation::NW);
            else return pare_neighbor->GetChild(Orientation::NE);
        }
    }
    else if(Direction::East == direction){
        if(node == parent->GetChild(Orientation::NW)) return parent->GetChild(Orientation::NE);
        if(node == parent->GetChild(Orientation::SW)) return parent->GetChild(Orientation::SE);

        QuadNode * pare_neighbor = Neighbor(parent, root, Direction::East);
        if(nullptr == pare_neighbor || !pare_neighbor->hasChild()) return pare_neighbor;
        else{
            if(node == parent->GetChild(Orientation::NE)) return pare_neighbor->GetChild(Orientation::NW);
            else return pare_neighbor->GetChild(Orientation::SW);
        }
    }
    else if(Direction::North == direction){
        if(node == parent->GetChild(Orientation::SW)) return parent->GetChild(Orientation::NW);
        if(node == parent->GetChild(Orientation::SE)) return parent->GetChild(Orientation::NE);

        QuadNode * pare_neighbor = Neighbor(parent, root, Direction::North);
        if(nullptr == pare_neighbor || !pare_neighbor->hasChild()) return pare_neighbor;
        else{
            if(node == parent->GetChild(Orientation::NW)) return pare_neighbor->GetChild(Orientation::SW);
            else return pare_neighbor->GetChild(Orientation::SE);
        }
    }
    else{
        if(node == parent->GetChild(Orientation::NE)) return parent->GetChild(Orientation::NW);
        if(node == parent->GetChild(Orientation::SE)) return parent->GetChild(Orientation::SW);

        QuadNode * pare_neighbor = Neighbor(parent, root, Direction::West);
        if(nullptr == pare_neighbor || !pare_neighbor->hasChild()) return pare_neighbor;
        else{
            if(node == parent->GetChild(Orientation::NW)) return pare_neighbor->GetChild(Orientation::NE);
            else return pare_neighbor->GetChild(Orientation::SE);
        }
    }
}

template <typename num_type, typename object, typename extent>
template <typename Condition>
inline void QuadTree<num_type, object, extent>::CreateSubNodesIf(QuadNode * node, Condition && condition)
{
    if(!node->hasChild() && condition(node)) node->CreateSubNodes();
    if(node->hasChild()){
        for(auto & child : node->GetChildren())
            CreateSubNodesIf(child, std::forward<Condition>(condition));
    }
}

template <typename num_type, typename object, typename extent>
template <typename T>
inline bool QuadTree<num_type, object, extent>::CreateSubBoxCondition() const
{
    return CreateSubBoxCondition(typename num_traits_float_or_int<T>::tag());
}

template <typename num_type, typename object, typename extent>
inline bool QuadTree<num_type, object, extent>::CreateSubBoxCondition(const num_floating_tag) const
{
    num_type epsilon = std::numeric_limits<num_type>::epsilon();
    if(m_bbox.Length() > 0 && m_bbox.Area() > epsilon) return true;
    return false;
}

template <typename num_type, typename object, typename extent>
inline bool QuadTree<num_type, object, extent>::CreateSubBoxCondition(const num_integer_tag) const
{
    if(m_bbox.Length() > 1 && m_bbox.Width() > 1) return true;
    return false;
}

}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_QUADTREE_HPP
