/**
 * @file Builder.hpp
 * @author bwu
 * @brief General tree building helper classes
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_TREE_BUILDER_HPP
#define GENERIC_TREE_BUILDER_HPP
#include <stack>
namespace generic{
///@brief tree structures and algorithms
namespace tree{

///@brief intermidiate store class for partition build
struct WorkItem
{
    size_t nodeIndex;
    size_t begin;
    size_t end;
    size_t depth;
    WorkItem() = default;
    WorkItem(size_t _nodeIndex, size_t _begin, size_t _end, size_t _depth)
        : nodeIndex(_nodeIndex), begin(_begin), end(_end), depth(_depth)
    {}

    size_t WorkSize() const { return end - begin; }
};

///@brief a single thread task spawner when build tree from top to down
struct TopDownTaskSpawner
{
    size_t taskSpawnThreshold = 1024;
    template <typename BuildTask>
    void RunTask(BuildTask & task, WorkItem item)
    {
        std::stack<WorkItem > stack;
        stack.emplace(item);
        while(!stack.empty())
        {
            WorkItem workItem = stack.top();
            stack.pop();

            auto moreItems = task.Build(workItem);
            if(moreItems.size()){
                auto iter = moreItems.begin();
                auto firstItem = *iter;
                iter++;
                for(; iter != moreItems.end(); ++iter){
                    if(firstItem.WorkSize() > iter->WorkSize())
                        std::swap(*iter, firstItem);
                    stack.push(*iter);
                }

                if(firstItem.WorkSize() > taskSpawnThreshold){
                    BuildTask subTask(task);
                    RunTask(subTask, firstItem);
                }
                else{
                    stack.push(firstItem);
                }
            }
        }
    }
};

}//namespace tree
}//namespace generic
#endif//GENERIC_TREE_BUILDER_HPP