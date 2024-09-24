/**
 * @file ThreadPool.hpp
 * @author bwu
 * @brief A header only thread pool implementation
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "Utility.hpp"
#include <boost/lockfree/queue.hpp>
#include <functional>
#include <future>
#include <atomic>
#include <thread>
#include <vector>
namespace generic{
///@brief thread related classes and functions
namespace thread {

class ThreadJoiner
{
public:
    explicit ThreadJoiner(std::vector<std::unique_ptr<std::thread> > & threads) : m_threads(threads){}
    ~ThreadJoiner()
    {
        Join();
    }

    void Join()
    {
        for(size_t i = 0; i < m_threads.size(); ++i){
            if(m_threads[i]->joinable())
                m_threads[i]->join();
        }
    }

private:
    std::vector<std::unique_ptr<std::thread> > & m_threads;
};

///@brief a thread pool class that holds and runs tasks parallelly
class ThreadPool
{
    using Task = std::function<void()>;
    using PoolQueue = boost::lockfree::queue<Task *, boost::lockfree::fixed_sized<false> >;
public:
    ///@brief constructs a pool with hardware supported threads
    ThreadPool();
    ///@brief constructs a poll with min(input threads, hardware supported threads) threads
    ThreadPool(size_t threads, size_t queueSize = 128);
    ThreadPool(ThreadPool & other) = delete;
    ThreadPool(const ThreadPool & other) = delete;
    ~ThreadPool();

    ThreadPool & operator= (const ThreadPool & pool) = delete;

    ///@brief get thread number used in this pool
    size_t Threads() const { return m_threads.size(); }
    ///@brief get current idel thread number in this pool
    size_t IdleThreads() const { return m_nWait.load(); }
    ///@brief check whether current pool is empty with task
    bool isEmpty() const { return m_queue.empty(); }
    ///@brief accesses thread by index
    std::thread & GetThread(size_t i) { return *(m_threads[i]); }
    ///@brief resizes thread number of this pool
    void Resize(size_t tryThreads);
    /**
     * @brief submit a task to the pool
     * 
     * @tparam FunctionType could be one of std::bind, std::function or lambda expression  
     * @param f Task function object
     * @return std::future of the FunctionType result 
     */
    template <typename FunctionType>
    std::future<std::invoke_result_t<FunctionType> >
    Submit(FunctionType && f);

    ///@brief pop one task from task queue
    Task PopTask();
    ///@brief wait for all tasks done
    void Wait();
    ///@brief wait for runing tasks done then stop
    void Stop();
    ///@brief empties the task queue
    void EmptyQueue();
    ///@brief gets current avaliable hardware thread number
    static size_t AvailableThreads(size_t tryThreads);

private:
    void Init_(size_t threads);
    void SetThread_(size_t thread);

private:
    std::vector<std::unique_ptr<std::thread> > m_threads;
    std::vector<std::shared_ptr<std::atomic_bool> > m_flags;
    mutable PoolQueue m_queue;
    std::atomic_bool m_bDone;
    std::atomic_bool m_bStop;
    std::atomic<size_t> m_nWait;
    mutable std::mutex m_mutex;
    std::condition_variable m_cond;
    ThreadJoiner m_joiner;
};

inline ThreadPool::ThreadPool()
 : m_queue(128), m_joiner(m_threads)
{
    Init_(std::thread::hardware_concurrency());
 }

inline ThreadPool::ThreadPool(size_t threads, size_t queueSize)
 : m_queue(queueSize), m_joiner(m_threads)
{
    Init_(threads);
}

inline ThreadPool::~ThreadPool()
{
    Wait();
}

inline void ThreadPool::Init_(size_t threads)
{
    m_nWait = 0; m_bStop = false; m_bDone = false;
    Resize(threads);
}

inline void ThreadPool::Resize(size_t tryThreads)
{
    if(tryThreads == 0) tryThreads = std::numeric_limits<size_t>::max();
    size_t nThreads = AvailableThreads(tryThreads);
    size_t currThreads = m_threads.size();
    if(!m_bStop && !m_bDone){
        if(currThreads <= nThreads){
            m_threads.resize(nThreads);
            m_flags.resize(nThreads);
        }

        for(size_t i = currThreads; i < nThreads; ++i){
            m_flags[i] = std::make_shared<std::atomic_bool>(false);
            SetThread_(i);
        }
    }
    else{
        for(size_t i = currThreads - 1; i >= nThreads; --i){
            *(m_flags[i]) = true;
            m_threads[i]->detach();
        }
        {
            std::unique_lock<std::mutex> lock(m_mutex);
            m_cond.notify_all();
        }
        m_threads.resize(nThreads);
        m_flags.resize(nThreads);
    }
}

inline void ThreadPool::EmptyQueue()
{
    Task * task;
    while(m_queue.pop(task))
        delete task;
}

inline typename ThreadPool::Task
ThreadPool::PopTask()
{
    Task * t = nullptr;
    m_queue.pop(t);
    std::unique_ptr<Task> task(t);

    Task rt;
    if(t) rt = *t;
    return rt;
}

inline void ThreadPool::SetThread_(size_t thread)
{
    std::shared_ptr<std::atomic_bool> flag(m_flags[thread]);
    auto f = [this, flag](){
        std::atomic_bool & bflag = *flag;
        Task * task;
        bool bPop = m_queue.pop(task);
        while(true){
            while(bPop){
                std::unique_ptr<Task> t(task);
                (*task)();

                if(bflag) return;
                else bPop = m_queue.pop(task);
            }
            
            std::unique_lock<std::mutex> lock(m_mutex);
            ++m_nWait;
            m_cond.wait(lock, [this, &task, &bPop, &bflag]()
            { bPop = m_queue.pop(task); return (bPop || m_bDone || bflag);});
            --m_nWait;

            if(!bPop) return;
        }
    };
    m_threads[thread].reset(new std::thread(f));
}

template <typename FunctionType>
inline std::future<std::invoke_result_t<FunctionType> >
ThreadPool::Submit(FunctionType && f)
{
    using result_type = std::invoke_result_t<FunctionType>;
    auto pkg = std::make_shared<std::packaged_task<result_type()> >(std::forward<FunctionType>(f));
    auto task = new std::function<void()>([pkg]{(*pkg)();});
    while(not m_queue.push(task));

    std::unique_lock<std::mutex> lock(m_mutex);
    m_cond.notify_one();
    
    return pkg->get_future();
}

inline void ThreadPool::Wait()
{
    if(m_bDone || m_bStop) return;
    m_bDone = true;
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cond.notify_all();
    }
    m_joiner.Join();
    EmptyQueue();
    m_threads.clear();
    m_flags.clear();
}

inline void ThreadPool::Stop()
{
    if(m_bStop) return;
    m_bStop = true;
    for(size_t i = 0; i < m_flags.size(); ++i)
        *(m_flags[i]) = true;
    EmptyQueue();
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_cond.notify_all();
    }
    m_joiner.Join();
    m_threads.clear();
    m_flags.clear();
}

inline size_t ThreadPool::AvailableThreads(size_t tryThreads)
{
    size_t currThreads = std::thread::hardware_concurrency();
    if(tryThreads > currThreads) tryThreads = currThreads;
    if(tryThreads == 0) tryThreads = 2;
    return tryThreads;
}
}//namespace thread
}//namespace generic