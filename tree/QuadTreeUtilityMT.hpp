#ifndef GENERIC_TREE_QUADTREEUTILITYMT_HPP
#define GENERIC_TREE_QUADTREEUTILITYMT_HPP
#include "QuadTree.hpp"
#include "thread/ThreadPool.hpp"
namespace generic{
namespace tree
{
template <typename object, typename Tree>
class QuadTreeBuilderMT
{ 
    using Pool = thread::ThreadPool;
    using Node = typename Tree::QuadNode;
    using Children = typename Tree::QuadChildren;
public:
    QuadTreeBuilderMT(Tree & tree, size_t threads = std::numeric_limits<size_t>::max())
     : m_tree(tree), m_thread(threads) {}

    void Build(std::list<object * > objs, size_t max_objs_to_build_sub);
    void BuildNode(Node * node, size_t max_objs_to_build_sub);
    void BuildNodeInThread(Node * node, size_t max_objs_to_build_sub);

private:
    Tree & m_tree;
    Pool m_thread;
};

template <typename object, typename Tree>
void QuadTreeBuilderMT<object, Tree>::Build(std::list<object* > objs, size_t max_objs_to_build_sub)
{
    m_tree.Build(std::move(objs));
    BuildNode(&m_tree, max_objs_to_build_sub);
    m_thread.Wait();
}

template <typename object, typename Tree>
void QuadTreeBuilderMT<object, Tree>::BuildNode(Node * node, size_t max_objs_to_build_sub)
{
    if(node->CreateSubNodes(max_objs_to_build_sub)){
        Children & children = node->GetChildren();
        for(auto * child : children){
            if(m_thread.isEmpty())
                BuildNodeInThread(child, max_objs_to_build_sub);
            else BuildNode(child, max_objs_to_build_sub);
        }
    }
}

template <typename object, typename Tree>
void QuadTreeBuilderMT<object, Tree>::BuildNodeInThread(Node * node, size_t max_objs_to_build_sub)
{
    m_thread.Submit(std::bind(&QuadTreeBuilderMT<object, Tree>::BuildNode, this, node, max_objs_to_build_sub));
}    
}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_QUADTREEUTILITYMT_HPP