/**
 * @file BuilderMT.hpp
 * @author bwu
 * @brief General tree building helper classes in multi-threads
 * @version 0.1
 * @date 2022-02-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once
#include "generic/thread/ThreadPool.hpp"
#include "Builder.hpp"
namespace generic{
namespace tree {

using generic::thread::ThreadPool;
///@brief a multi threads task spawner when build tree from top to down
struct TopDownTaskSpawnerMT
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
                    if(m_thread.isEmpty())
                        RunTaskInThread(subTask, firstItem);
                    else RunTask(subTask, firstItem);
                }
                else{
                    stack.push(firstItem);
                }
            }
        }
    }
    
    template <typename BuildTask>
    void RunTaskInThread(BuildTask & task, WorkItem item)
    {
        m_thread.Submit(std::bind(&TopDownTaskSpawnerMT::RunTask<BuildTask>, this, task, item));
    }

private:
    ThreadPool m_thread;
};
}//namespace tree
}//namespace generic
