/**
 * @file TaskFlow.hpp
 * @author bwu
 * @brief A header only task flow library that run tasks based on dependency parallelly
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_THREAD_TASKFLOW_HPP
#define GENERIC_THREAD_TASKFLOW_HPP
#include "ThreadPool.hpp"
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graphviz.hpp>
#include <condition_variable>
#include <atomic>
#include <vector>
#include <deque>
#include <list>
namespace generic {
namespace thread {
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
    ///@brief constructs a task node with node id and funtion object to be executed
    template<typename FunctionType>
    TaskNode(size_t id, FunctionType && func);

    ///@brief add precede dependency for this node, input tasks will be executed after this task
    template <typename... Args>
    void Precede(TaskNode * task, Args &&... args);

    ///@brief add success dependency for this node, this task will be executed after input tasks
    template <typename... Args>
    void Success(TaskNode * task, Args &&... args);

    void Precede(){}
    void Success(){}

private:
    size_t m_id;
    std::list<TaskNode * > m_successors;
    std::list<TaskNode * > m_dependents;
    std::unique_ptr<FunctionWrapper > m_work;
};

template<typename FunctionType>
inline TaskNode::TaskNode(size_t id, FunctionType && func)
 : m_id(id), m_work(new FunctionWrapper(std::move(func))) {}

template <typename... Args>
inline void TaskNode::Precede(TaskNode * task, Args && ...args)
{
    m_successors.push_back(task);
    Precede(std::forward<Args &&>(args)...);
}

template <typename... Args>
inline void TaskNode::Success(TaskNode * task, Args && ...args)
{
    m_dependents.push_back(task);
    Success(std::forward<Args &&>(args)...);
}

///@brief task flow class for registering and managing tasks
class TaskFlow
{
    using TaskGraph = boost::adjacency_list<
                      boost::vecS,
                      boost::vecS,
                      boost::bidirectionalS,
                      boost::property<boost::vertex_index_t, size_t >
                      >;
    using Vertex = boost::graph_traits<TaskGraph>::vertex_descriptor;

    class TaskLabelWriter
    {
    public:
        TaskLabelWriter(const std::unordered_map<size_t, std::string> & labels) : m_labels(labels){}

        void operator() (std::ostream & out, const Vertex & v) const
        {
            if(m_labels.count(v)){ out << "[label=\"" << m_labels.at(v) << "\"]"; }
            else { out << "[label=\"" << "task" << std::to_string(v) << "\"]"; }
        }
    private:
        const std::unordered_map<size_t, std::string> & m_labels;
    };

    friend Executor;
    friend Dispatcher;
public:
    TaskFlow(){}
    ~TaskFlow(){}
    TaskFlow(TaskFlow & other) = delete;
    TaskFlow(const TaskFlow & other) = delete;
    TaskFlow & operator= (TaskFlow & other) = delete;

    ///@brief accesses task node by index
    TaskNode& operator[] (size_t i) {return *m_tasks[i];}

    /**
     * @brief submit a task and get the internal task node
     * @tparam FunctionType could be one of std::bind, std::function or lambda expression
     * @param f Task function object
     * @param label task label name
     * @return a pair with internal task node and std::future of the FunctionType result
     */
    template <typename FunctionType>
    std::pair<TaskNode *,
    std::future<typename std::result_of<FunctionType()>::type> >
    Submit(FunctionType && f, std::string label = std::string());

    /**
     * @brief emplace a task and get the internal task node
     * 
     * @tparam FunctionType could be one of std::bind, std::function or lambda expression 
     * @param f Task function object
     * @param label task label name
     * @return TaskNode* internal task node
     */
    template <typename FunctionType>
    TaskNode * Emplace(FunctionType && f, std::string label = std::string());

    ///@brief build a direct task flow graph based on added tasks dependency
    void BuildTaskGraph();
    ///@brief print *.dot format task graph to out stream
    void PrintTaskGraph(std::ostream & out);
private:
    void TopologicalSort();
    std::list<TaskNode * > Successors(TaskNode & node);
    std::list<TaskNode * > Dependents(TaskNode & node);

private:
    std::deque<Vertex> m_sequence;
    std::unique_ptr<TaskGraph> m_graph;
    std::vector<std::unique_ptr<TaskNode> > m_tasks;
    std::unordered_map<size_t, std::string> m_taskLabels;
    std::vector<std::unique_ptr<std::atomic<size_t> > > m_dependents;
};

template <typename FunctionType>
std::pair<TaskNode *,
std::future<typename std::result_of<FunctionType()>::type> >
inline TaskFlow::Submit(FunctionType && f, std::string label)
{
    using result_type = typename std::result_of<FunctionType()>::type;
    std::packaged_task<result_type()> task(std::move(f));
    std::future<result_type> res(task.get_future());

    size_t id = m_tasks.size();
    m_tasks.emplace_back(new TaskNode(id, std::move(task)));

    if(label.empty()) label = "task" + std::to_string(id);
    m_taskLabels.insert(std::make_pair(id, label));
    return std::make_pair(m_tasks.back().get(), std::move(res));
}

template <typename FunctionType>
inline TaskNode *  TaskFlow::Emplace(FunctionType && f, std::string label)
{
    size_t id = m_tasks.size();
    m_tasks.emplace_back(new TaskNode(id, std::move(f)));

    if(label.empty()) label = "task" + std::to_string(id);
    m_taskLabels.insert(std::make_pair(id, label));
    return m_tasks.back().get();
}

inline std::list<TaskNode * > TaskFlow::Successors(TaskNode & node)
{
    if(nullptr == m_graph) BuildTaskGraph();

    std::list<TaskNode * > successors;
    auto ee = boost::out_edges(node.m_id, *m_graph);
    for(auto e = ee.first; e != ee.second; ++e){
        Vertex v = boost::target(*e, *m_graph);
        successors.push_back(&(*m_tasks[v]));
    }
    return successors;
}

inline std::list<TaskNode * > TaskFlow::Dependents(TaskNode & node)
{
    if(nullptr == m_graph) BuildTaskGraph();

    std::list<TaskNode * > dependents;
    auto ee = boost::in_edges(node.m_id, *m_graph);
    for(auto e = ee.first; e != ee.second; ++e){
        Vertex v = boost::source(*e, *m_graph);
        dependents.push_back(&(*m_tasks[v]));
    }
    return dependents;
}

inline void TaskFlow::BuildTaskGraph()
{
    size_t taskSize = m_tasks.size();
    m_graph.reset(new TaskGraph(taskSize));
    for(size_t i = 0; i < taskSize; ++i){
        auto & node = m_tasks[i];
        for(TaskNode * successor : node->m_successors)
            boost::add_edge(i, successor->m_id, *m_graph);
        for(TaskNode * dependent : node->m_dependents)
            boost::add_edge(dependent->m_id, i, *m_graph);
    }

    m_dependents.resize(taskSize);
    for(size_t i = 0; i < taskSize; ++i)
        m_dependents[i].reset(new std::atomic<size_t>(boost::in_degree(i, *m_graph)));
}

inline void TaskFlow::PrintTaskGraph(std::ostream & out)
{
    if(nullptr == m_graph) BuildTaskGraph();
    boost::write_graphviz(out, *m_graph, TaskLabelWriter(m_taskLabels));
}

inline void TaskFlow::TopologicalSort()
{
    m_sequence.clear();
    boost::topological_sort(*m_graph, std::front_inserter(m_sequence));
}

///@brief dispatcher class that dispatch tasks based on dependency
class Dispatcher
{
public:
    explicit Dispatcher(TaskFlow & flow) : m_flow(flow){}
    ~Dispatcher(){}
    
    void Dispatch(size_t threads);
private:
    void DispatchOne(TaskNode & node);
private:
    TaskFlow & m_flow;
    mutable std::mutex m_mutex;
    std::condition_variable m_cond;
};

inline void Dispatcher::Dispatch(size_t threads)
{
    ThreadPool pool(threads);
    while(!m_flow.m_sequence.empty()){
        TaskFlow::Vertex v = m_flow.m_sequence.front();
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cond.wait(lock, [&]{ return m_flow.m_dependents[v]->load() == 0; });
        
        size_t count = 0;
        auto iter = m_flow.m_sequence.begin();
        for(; iter != m_flow.m_sequence.end() && count < threads;){
            if(m_flow.m_dependents[*iter]->load() == 0){
                TaskNode & node = m_flow[*iter];
                pool.Submit(std::bind(&Dispatcher::DispatchOne, this, std::ref(node)));
                iter = m_flow.m_sequence.erase(iter);
                count++;
            }
            else{ ++iter; }
        }
    }
}

inline void Dispatcher::DispatchOne(TaskNode & node)
{
    (*node.m_work)();
    std::list<TaskNode * > successors = m_flow.Successors(node);
    for(auto * successor : successors){
        m_flow.m_dependents[successor->m_id]->fetch_sub(1);
    }
    m_cond.notify_all();
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
     * @return false if any task has circular dependency 
     */
    bool Run(TaskFlow & flow);
private:
    size_t m_threads;
};

inline bool Executor::Run(TaskFlow & flow)
{
    flow.BuildTaskGraph();
    try { flow.TopologicalSort(); }
    catch (...) { return false; }

    Dispatcher dispatcher(flow);
    dispatcher.Dispatch(m_threads);
    return true;
}

}//namespace taskflow
}//namespace thread
}//namespace generic
#endif//GENERIC_THERAD_TASKFLOW_HPP