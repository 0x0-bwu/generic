#ifndef GENERIC_TREE_RECTTREE_HPP
#define GENERIC_TREE_RECTTREE_HPP
#include "geometry/Box.hpp"
#include <boost/geometry.hpp>
#include <boost/geometry/core/point_type.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <iostream>
#include <utility>
#include <list>
namespace generic{
namespace tree {
//this is a boost R(ectangle) tree wapper for the convenient usage
using generic::geometry::Box2D;
using generic::geometry::Point2D;
namespace bg = boost::geometry;
namespace bgid = bg::index::detail;

template <typename num_type, typename object>
class RectNode
{
public:
    typedef std::list<RectNode<num_type, object> *> RectChildren;
    int m_level = 0;
    Box2D<num_type > m_bbox;
    std::list<object * > m_objs;
    RectNode * m_parent = nullptr;
    RectChildren m_children;

    virtual ~RectNode(){ Clear(); }
    void Build(std::list<object* > && objs);

    int GetLevel() const { return m_level; }
    void SetBBox(const Box2D<num_type > & bbox) { m_bbox = bbox; }
    const Box2D<num_type > & GetBBox() const { return m_bbox; }

    std::list<object * > & GetObjs() { return m_objs; }
    const std::list<object * > & GetObjs() const { return m_objs; }

    bool isRoot() const { return nullptr == m_parent; }
    RectNode * GetParent() { return m_parent; }
    bool hasChild() const { return m_children.size() != 0 ; }
    RectChildren & GetChildren() { return m_children; }
    const RectChildren & GetChildren() const { return m_children; }
    void GetAllObjects(std::list<object * > & allObjs) const;
    void GetAllObjects(std::list<const object * > & allObjs) const;
    static void GetAllNodes(RectNode<num_type, object>  * node, std::list<RectNode<num_type, object> * > & allNodes);
    static void GetAllNodesByLevel(RectNode<num_type, object>  * node, std::map<int, std::list<RectNode<num_type, object> * > > & allNodes);
    static void GetTopNodesAreaLessThan(num_type area, RectNode<num_type, object> * node, std::list<RectNode<num_type, object> * > & allNodes);

    void Clear();
};

template <typename num_type, typename object>
void RectNode<num_type, object>::Build(std::list<object * > && objs)
{
    Clear();
    m_objs = objs;
}

template <typename num_type, typename object>
void RectNode<num_type, object>::GetAllObjects(std::list<object * > & allObjs) const
{
    typename std::list<object* >::const_iterator iter1 = m_objs.begin();
    for(; iter1 != m_objs.end(); ++iter1) allObjs.push_back(*iter1);

    if(hasChild()){
        auto iter = m_children.begin();
        for(; iter != m_children.end(); ++iter)
            (*iter)->GetAllObjects(allObjs);
    }
}

template <typename num_type, typename object>
void RectNode<num_type, object>::GetAllObjects(std::list<const object * > & allObjs) const
{
    typename std::list<object * >::const_iterator iter1 = m_objs.begin();
    for(; iter1 != m_objs.end(); ++iter1) allObjs.push_back(*iter1);

    if(hasChild()){
        auto iter = m_children.begin();
        for(; iter != m_children.end(); ++iter)
            (*iter)->GetAllObjects(allObjs);
    }
}

template <typename num_type, typename object>
void RectNode<num_type, object>::GetAllNodes(RectNode<num_type, object>  * node, std::list<RectNode<num_type, object> * > & allNodes)
{
    allNodes.push_back(node);
    if(!node->hasChild()) return;
    const auto & children = node->GetChildren();
    for(auto & child : children) GetAllNodes(child, allNodes);
}

template <typename num_type, typename object>
void RectNode<num_type, object>::GetAllNodesByLevel(RectNode<num_type, object> * node, std::map<int, std::list<RectNode<num_type, object> * > > & allNodes)
{
    int topLevel = node->GetLevel();
    allNodes[topLevel].push_back(node);
    if(!node->hasChild()) return;
    const auto & children = node->GetChildren();
    for(auto & child : children) GetAllNodesByLevel(child, allNodes);
}

template <typename num_type, typename object>
void RectNode<num_type, object>::GetTopNodesAreaLessThan(num_type area, RectNode<num_type, object> * node, std::list<RectNode<num_type, object> * > & allNodes)
{
    if(node->m_bbox.Area() >= area){
        if(!node->hasChild()) return;
        const auto & children = node->GetChildren();
        for(auto & child : children) GetTopNodesAreaLessThan(area, child, allNodes);    
    }
    else allNodes.push_back(node);
}

template <typename num_type, typename object>
void RectNode<num_type, object>::Clear()
{
    auto iter = m_children.begin();
    for(; iter != m_children.end(); ++iter){
        if(*iter) delete *iter;
    }
    m_children.clear();
    m_objs.clear();
}

template <typename num_type, typename object, typename extent>
class RectTree : public RectNode<num_type, object>
{
    typedef bg::model::d2::point_xy<num_type, bg::cs::cartesian> RtPoint;
    typedef bg::model::box<RtPoint> RtBox;
    typedef std::pair<RtBox, object * > RtValue;
    typedef bg::index::rtree<RtValue, bg::index::dynamic_rstar> RtreeProxy;
    typedef bgid::rtree::utilities::view<RtreeProxy> RtProxyView;

    template <typename Value, typename Options, typename Translator, typename Box, typename Allocators>
    class RtProxyVisitor
            : public bgid::rtree::visitor<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag, true>::type
    {
        typedef typename bgid::rtree::internal_node<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type RtProxyNode;
        typedef typename bgid::rtree::leaf<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type RtProxyLeaf;

    public:
        RtProxyVisitor(RectNode<num_type, object> * root) : m_root(root) { m_currNode = m_root; }

        void operator()(const RtProxyNode & proxyNode)
        {
            typedef typename bgid::rtree::elements_type<RtProxyNode>::type RtProxyNodeEle;
            const RtProxyNodeEle & elements = bgid::rtree::elements(proxyNode);

            RectNode<num_type, object> * pareNode = m_currNode;
            auto iter = elements.begin();
            for(; iter != elements.end(); ++iter){
                RectNode<num_type, object> * node = new RectNode<num_type, object>;
                HandleBox(node, pareNode, iter->first);
                m_currNode = node;
                bgid::rtree::apply_visitor(*this, *(iter->second));
            }
        }

        void operator()(const RtProxyLeaf & leaf)
        {
            typedef typename bgid::rtree::elements_type<RtProxyLeaf>::type RtProxyLeafEle;
            const RtProxyLeafEle & elements = bgid::rtree::elements(leaf);

            auto iter = elements.begin();
            for(; iter != elements.end(); ++iter){
                HandleValue(m_currNode, iter->second);
            }
        }

        void HandleBox(RectNode<num_type, object> * node, RectNode<num_type, object> * parent, const Box & box)
        {
            node->m_level = parent->m_level + 1;
            node->m_parent = parent;
            parent->m_children.push_back(node);
            node->m_bbox[0][0] = bg::get<bg::min_corner, 0>(box);
            node->m_bbox[0][1] = bg::get<bg::min_corner, 1>(box);
            node->m_bbox[1][0] = bg::get<bg::max_corner, 0>(box);
            node->m_bbox[1][1] = bg::get<bg::max_corner, 1>(box);
        }

        void HandleValue(RectNode<num_type, object> * node, object * obj)
        {
            node->m_objs.push_back(obj);
        }

        RectNode<num_type, object> * m_root;
        RectNode<num_type, object> * m_currNode;
    };

public:
    RectTree() : m_proxy(nullptr) {}
    ~RectTree()
    {
        if(m_proxy) delete m_proxy;
        m_proxy = nullptr;
        RectNode<num_type, object>::Clear();
    }

    void Build(std::list<object* > objs, size_t max_objs_to_build_sub);
    void GetTopNodesAreaLessThan(num_type area, std::list<RectNode<num_type, object> * > & nodes);

private:
    void BuildProxy_(const std::list<object* > & objs, size_t max_objs_to_build_sub);
    void BuildSelf_();

    RtreeProxy * m_proxy;
};

template <typename num_type, typename object, typename extent>
void RectTree<num_type, object, extent>::Build(std::list<object* > objs, size_t max_objs_to_build_sub)
{
    RectNode<num_type, object>::Clear();

    if(0 == max_objs_to_build_sub){
        RectNode<num_type, object>::Build(std::move(objs));
        return;
    }

    BuildProxy_(objs, max_objs_to_build_sub);
    BuildSelf_();
}

template <typename num_type, typename object, typename extent>
void RectTree<num_type, object, extent>::BuildProxy_(const std::list<object* > & objs, size_t max_objs_to_build_sub)
{
    extent ext;
    if(m_proxy) delete m_proxy;
    bg::index::dynamic_rstar para(max_objs_to_build_sub);
    m_proxy = new RtreeProxy(para);

    for(object * obj : objs){
        Box2D<num_type> bbox = ext(*obj);
        RtBox rtBox(RtPoint(bbox[0][0], bbox[0][1]), RtPoint(bbox[1][0], bbox[1][1]));
        m_proxy->insert(std::make_pair(rtBox, obj));
    }
}

template <typename num_type, typename object, typename extent>
void RectTree<num_type, object, extent>::BuildSelf_()
{
    RtProxyView rtv(*m_proxy);
    RtProxyVisitor<
            typename RtProxyView::value_type,
            typename RtProxyView::options_type,
            typename RtProxyView::translator_type,
            typename RtProxyView::box_type,
            typename RtProxyView::allocators_type
            > visitor(this);
    rtv.apply_visitor(visitor);
}

template <typename num_type, typename object, typename extent>
void RectTree<num_type, object, extent>::GetTopNodesAreaLessThan(num_type area, std::list<RectNode<num_type, object> * > & nodes)
{
    RectNode<num_type, object>::GetTopNodesAreaLessThan(area, this, nodes);
}

}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_RECTTREE_HPP