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

class TaskNode
{
    friend TaskFlow;
    friend Executor;
    friend Dispatcher;
public:
    template<typename FunctionType>
    TaskNode(size_t id, FunctionType && func);

    template <typename... Args>
    void Precede(TaskNode * task, Args &&... args);
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

    TaskNode& operator[] (size_t i) {return *m_tasks[i];}

    template <typename FunctionType>
    std::pair<TaskNode *,
    std::future<typename std::result_of<FunctionType()>::type> >
    Submit(FunctionType && f, std::string label = std::string());

    template <typename FunctionType>
    TaskNode * Emplace(FunctionType && f, std::string label = std::string());

    void BuildTaskGraph();
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

class Executor
{
public:
    Executor(size_t threads) : m_threads(threads) {}
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