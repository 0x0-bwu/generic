/**
 * @file TaskFlow.hpp
 * @author bwu
 * @brief A header only task flow library that run tasks based on dependency parallelly
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "ThreadPool.hpp"
#include <condition_variable>
#include <unordered_map>
#include <unordered_set>
#include <optional>

namespace generic {
namespace thread {
///@brief namespace of task flow library
namespace taskflow {

class TaskFlow;
class Executor;
class Dispatcher;
class TaskLabelWriter;

///@brief task node that hold the task to be executed and dependency relationships with other task nodes
class TaskNode
{
    friend TaskFlow;
    friend Executor;
    friend Dispatcher;
public:
    ///@brief constructs a task node with funtion object and label(optional) to be executed
    template<typename FunctionType>
    TaskNode(FunctionType && func, std::optional<std::string> label);

    ///@brief add precede dependency for this node, input tasks will be executed after this task
    template <typename... Args>
    void Precede(TaskNode * task, Args &&... args);

    ///@brief add success dependency for this node, this task will be executed after input tasks
    template <typename... Args>
    void Success(TaskNode * task, Args &&... args);

/**
 * @brief Brief description of Precede.
 * @return void
 */
    void Precede(){}
/**
 * @brief Brief description of Success.
 * @return void
 */
    void Success(){}

private:
    void PrecedeImp(TaskNode * task);
    void SuccessImp(TaskNode * task);

private:
    std::unordered_set<TaskNode * > m_successors;
    std::unordered_set<TaskNode * > m_dependents;
    std::unique_ptr<FunctionWrapper > m_work;
    std::optional<std::string> m_label;
    std::atomic<size_t> m_depCount{0};
};

template<typename FunctionType>
inline TaskNode::TaskNode(FunctionType && func, std::optional<std::string> label)
 : m_work(new FunctionWrapper(std::move(func))), m_label(std::move(label)) {}

template <typename... Args>
inline void TaskNode::Precede(TaskNode * task, Args && ...args)
{
    this->PrecedeImp(task);
    task->SuccessImp(this);
    Precede(std::forward<Args &&>(args)...);
}

template <typename... Args>
inline void TaskNode::Success(TaskNode * task, Args && ...args)
{
    this->SuccessImp(task);
    task->PrecedeImp(this);
    Success(std::forward<Args &&>(args)...);
}

inline void TaskNode::PrecedeImp(TaskNode * task)
{
    m_successors.emplace(task);
}
    
inline void TaskNode::SuccessImp(TaskNode * task)
{
    if(m_dependents.emplace(task).second)
        m_depCount += 1;
}

///@brief task flow class for registering and managing tasks
class TaskFlow
{
    friend Executor;
    friend Dispatcher;
public:
    TaskFlow() = default;
    ~TaskFlow() = default;
    TaskFlow(TaskFlow & other) = delete;
    TaskFlow(const TaskFlow & other) = delete;
    TaskFlow & operator= (TaskFlow & other) = delete;

    /**
     * @brief submit a task and get the internal task node
     * @tparam FunctionType could be one of std::bind, std::function or lambda expression
     * @param f Task function object
     * @param label task label name
     * @return a pair with internal task node and std::shared_future of the FunctionType result
     */
    template <typename FunctionType>
    std::pair<TaskNode *,
    std::shared_future<std::invoke_result_t<FunctionType> > >
    Submit(FunctionType && f, std::optional<std::string> label = {});

    /**
     * @brief emplace a task and get the internal task node
     * 
     * @tparam FunctionType could be one of std::bind, std::function or lambda expression 
     * @param f Task function object
     * @param label task label name
     * @return TaskNode* internal task node
     */
    template <typename FunctionType>
    TaskNode * Emplace(FunctionType && f, std::optional<std::string> label = {});

    ///@brief print *.dot format task graph to out stream
    void PrintTaskGraph(std::ostream & out);

private:
    const std::unordered_set<TaskNode * > & Successors(TaskNode & node) const { return node.m_successors; }
    const std::unordered_set<TaskNode * > & Dependents(TaskNode & node) const { return node.m_dependents; }

private:
    bool m_executed = false;
    std::vector<std::unique_ptr<TaskNode> > m_tasks;
};

template <typename FunctionType>
std::pair<TaskNode *,
std::shared_future<std::invoke_result_t<FunctionType> > >
inline TaskFlow::Submit(FunctionType && f, std::optional<std::string> label)
{
    using ResultType = std::invoke_result_t<FunctionType>;
    std::packaged_task<ResultType()> task(std::move(f));
    std::shared_future<ResultType> res(task.get_future());

    m_tasks.emplace_back(new TaskNode(std::move(task), std::move(label)));
    return std::make_pair(m_tasks.back().get(), std::move(res));
}

template <typename FunctionType>
inline TaskNode *  TaskFlow::Emplace(FunctionType && f, std::optional<std::string> label)
{
    m_tasks.emplace_back(new TaskNode(std::move(f), std::move(label)));
    return m_tasks.back().get();
}

inline void TaskFlow::PrintTaskGraph(std::ostream & out)
{
    if(out.bad()) return;

    const char sp(' ');
    const char eol('\n');
    std::string indent(2, sp);
    out << "digraph g {" << eol;
    
    std::unordered_map<TaskNode *, size_t> idMap;
    for(size_t i = 0; i < m_tasks.size(); ++i) {
        idMap.emplace(m_tasks[i].get(), i);
        out << indent << 'n' << std::to_string(i) << "[label=\"" << m_tasks[i]->m_label.value_or("task" + std::to_string(i)) << "\"];" << eol;
    }

    for(size_t i = 0; i < m_tasks.size(); ++i) {
        for(auto * dep : m_tasks[i]->m_successors) {
            out << indent << 'n' << std::to_string(i) << sp << "->" << sp << 'n' << std::to_string(idMap.at(dep)) << ';' << eol;
        }
    }

    out << '}' << eol;
}

///@brief dispatcher class that dispatch tasks based on dependency
class Dispatcher
{
public:
/**
 * @brief Brief description of Dispatcher.
 * @param m_flow(flow
 * @return explicit
 */
    explicit Dispatcher(TaskFlow & flow) : m_flow(flow){}
    ~Dispatcher(){}
    
    void Dispatch(size_t threads);
private:
    void DispatchOne(TaskNode * node, ThreadPool * pool);
private:
    TaskFlow & m_flow;
    mutable std::mutex m_mutex;
    std::condition_variable m_cond;
    std::atomic<size_t> m_nTasks{0};
};

inline void Dispatcher::Dispatch(size_t threads)
{
    ThreadPool pool(threads);
    m_nTasks.store(m_flow.m_tasks.size());
    
    std::vector<TaskNode *> tasks;
    for(auto & task : m_flow.m_tasks) {
        if(task->m_depCount == 0)
            tasks.emplace_back(task.get());
    }

    for(auto * task : tasks)
        pool.Submit(std::bind(&Dispatcher::DispatchOne, this, task, &pool));
    
    std::unique_lock<std::mutex> lock(m_mutex);
    m_cond.wait(lock, [&]{ return m_nTasks.load() == 0; });
}

inline void Dispatcher::DispatchOne(TaskNode * node, ThreadPool * pool)
{
    (*node->m_work)();
    m_nTasks.fetch_sub(1);
    for(auto * successor : node->m_successors) {
        if(auto dep = successor->m_depCount.fetch_sub(1); dep == 1)
            pool->Submit(std::bind(&Dispatcher::DispatchOne, this, successor, pool));
    }
    m_cond.notify_one();
}

///@brief executor class to run taskflow in parallel
class Executor
{
public:
    ///@brief constructs a executor with specified thead numbers
    Executor(size_t threads) : m_threads(threads) {}

    /**
     * @brief execute input task flow `flow` parallelly in setting threads
     * @param[in] flow the take flow to be executed
     * @return true when task flow executed successfully
     * @return false if flow has been executed once or has circular dependency when the check enabled
     */
    bool Run(TaskFlow & flow);

    ///@brief set if check task cyclic dependency before execute
    void SetCheckCyclicDependency(bool check);

private:
    bool m_cycleCheck = true;
    size_t m_threads;
};

inline bool Executor::Run(TaskFlow & flow)
{
    if(flow.m_executed) return false;
    if(m_cycleCheck) {
        //wbtest, todo
    }

    flow.m_executed = true;
    Dispatcher dispatcher(flow);
    dispatcher.Dispatch(m_threads);
    return true;
}

inline void Executor::SetCheckCyclicDependency(bool check)
{
    m_cycleCheck = check;
} 

}//namespace taskflow
}//namespace thread
}//namespace generic
