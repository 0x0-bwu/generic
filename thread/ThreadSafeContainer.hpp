/**
 * @file ThreadSafeContainer.hpp
 * @author bwu
 * @brief Thread-safe container implementations
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include <mutex>
#include <queue>
#include <atomic>
#include <condition_variable>
namespace generic{
namespace thread {

///@brief a thread-safe task queue that can execute pop, push and other function simultaneously
template <typename Task>
class ThreadSafeQueue
{
    struct Node
    {
        std::shared_ptr<Task> task;
        std::unique_ptr<Node> next;
    };

public:
    ThreadSafeQueue();
    ThreadSafeQueue(const ThreadSafeQueue & queue) = delete;
    ThreadSafeQueue & operator=(const ThreadSafeQueue & queue) = delete;

    ///@brief pushes a task into the queue
    void Push(Task task);
    ///@brief tries to pop a task from the queue, returns false if queue is empty
    bool TryPop(Task & task);
    ///@brief waits for and pops a task from the queue
    void WaitAndPop(Task & task);
    ///@brief tries to pop a task from the queue, returns nullptr if queue is empty
    std::shared_ptr<Task> TryPop();
    ///@brief waits for and pops a task from the queue
    std::shared_ptr<Task> WaitAndPop();
    ///@brief checks if the queue is empty
    bool Empty();

private:
    Node * GetTail_();
    std::unique_ptr<Node> PopHead_();
    std::unique_lock<std::mutex> WatiForTask_();
    std::unique_ptr<Node> WaitPopHead_();
    std::unique_ptr<Node> WaitPopHead_(Task & task);
    std::unique_ptr<Node> TryPopHead_();
    std::unique_ptr<Node> TryPopHead_(Task & task);

private:
    mutable std::mutex m_headMutex;
    std::unique_ptr<Node> m_head;
    mutable std::mutex m_tailMutex;
    Node * m_tail = nullptr;
    std::condition_variable m_cond;
};

template <typename Task>
ThreadSafeQueue<Task>::ThreadSafeQueue() : m_head(new Node)
{
    m_tail = m_head.get();
}

template <typename Task>
inline void ThreadSafeQueue<Task>::Push(Task task)
{
    std::shared_ptr<Task> t(std::make_shared<Task>(std::move(task)));
    std::unique_ptr<Node> node(new Node);
    {
        std::lock_guard<std::mutex> lock(m_tailMutex);
        m_tail->task = t;
        Node * newTail = node.get();
        m_tail->next = std::move(node);
        m_tail = newTail;
    }
    m_cond.notify_one();
}

template <typename Task>
inline std::shared_ptr<Task> ThreadSafeQueue<Task>::TryPop()
{
    std::unique_ptr<Node> preHead = TryPopHead_();
    return preHead ? preHead->task : std::shared_ptr<Task>();
}

template <typename Task>
inline bool ThreadSafeQueue<Task>::TryPop(Task & task)
{
    const std::unique_ptr<Node> pre = TryPopHead_(task);
    return pre != nullptr;
}

template <typename Task>
inline std::shared_ptr<Task> ThreadSafeQueue<Task>::WaitAndPop()
{
    const std::unique_ptr<Node> pre = WaitPopHead_();
    return pre->task;
}

template <typename Task>
inline void ThreadSafeQueue<Task>::WaitAndPop(Task & task)
{
    WaitPopHead_(task);
}

template <typename Task>
inline bool ThreadSafeQueue<Task>::Empty()
{
    std::lock_guard<std::mutex> lock(m_headMutex);
    return (m_head.get() == GetTail_());
}

template <typename Task>
inline typename ThreadSafeQueue<Task>::Node * ThreadSafeQueue<Task>::GetTail_()
{
    std::lock_guard<std::mutex> lock(m_tailMutex);
    return m_tail;
}

template <typename Task>
inline std::unique_ptr<typename ThreadSafeQueue<Task>::Node > ThreadSafeQueue<Task>::PopHead_()
{
    std::unique_ptr<Node> preHead = std::move(m_head);
    m_head = std::move(preHead->next);
    return preHead;
}

template <typename Task>
inline std::unique_lock<std::mutex> ThreadSafeQueue<Task>::WatiForTask_()
{
    std::unique_lock<std::mutex> lock(m_headMutex);
    m_cond.wait(lock, [&]{ return m_head.get() != GetTail_(); });
    return std::move(lock);
}

template <typename Task>
inline std::unique_ptr<typename ThreadSafeQueue<Task>::Node > ThreadSafeQueue<Task>::WaitPopHead_()
{
    std::unique_lock<std::mutex> lock(WatiForTask_());
    return PopHead_();
}

template <typename Task>
inline std::unique_ptr<typename ThreadSafeQueue<Task>::Node > ThreadSafeQueue<Task>::WaitPopHead_(Task & task)
{
    std::unique_lock<std::mutex> lock(WatiForTask_());
    task = std::move(*m_head->task);
    return PopHead_();
}

template <typename Task>
inline std::unique_ptr<typename ThreadSafeQueue<Task>::Node > ThreadSafeQueue<Task>::TryPopHead_()
{
    std::lock_guard<std::mutex> lock(m_headMutex);
    if(m_head.get() == GetTail_()) return std::unique_ptr<Node>();
    return PopHead_();
}

template <typename Task>
inline std::unique_ptr<typename ThreadSafeQueue<Task>::Node > ThreadSafeQueue<Task>::TryPopHead_(Task & task)
{
    std::lock_guard<std::mutex> lock(m_headMutex);
    if(m_head.get() == GetTail_()) return std::unique_ptr<Node>();
    task = std::move(*m_head->task);
    return PopHead_();
}

template <typename Task>
class LockFreeQueue
{
    struct Node;
    struct CountedNodePtr
    {
        int externalCount = 0;
        Node * node = nullptr;
    };

    struct NodeCounter
    {
        unsigned internalCount:30;
        unsigned externalCounters:2;
    };

    struct Node
    {
        std::atomic<Task * > task;
        std::atomic<NodeCounter> count;
        CountedNodePtr next;

        Node()
        {
            NodeCounter c;
            c.internalCount = 0;
            c.externalCounters = 2;
            count.store(c);
        }

        void ReleaseRef()
        {
            NodeCounter preCounter = count.load(std::memory_order_relaxed);
            NodeCounter newCounter;
            do {
                newCounter = preCounter;
                --newCounter.internalCount;
            } while(!count.compare_exchange_strong(preCounter, newCounter, std::memory_order_acquire, std::memory_order_relaxed));
            if(!newCounter.internalCount && !newCounter.externalCounters) { delete this; }
        }
    };

public:
    LockFreeQueue();
    ~LockFreeQueue();
    LockFreeQueue(const LockFreeQueue & queue) = delete;
    LockFreeQueue & operator=(const LockFreeQueue & queue) = delete;

    std::unique_ptr<Task> Pop();
    void Push(Task task);

private:
    void SetNewTail_(CountedNodePtr & preTail, const CountedNodePtr & newTail);
    static void IncreaseExternalCount_(std::atomic<CountedNodePtr> & counter, CountedNodePtr & preCounter);
    static void FreeExternalCount_(CountedNodePtr & preNodePtr);

private:
    std::atomic<Node * > m_head;
    std::atomic<Node * > m_tail;
};

template <typename Task>
LockFreeQueue<Task>::LockFreeQueue() : m_head(new Node)
{
    m_tail.store(m_head.load());
}

template <typename Task>
LockFreeQueue<Task>::~LockFreeQueue()
{
    while(Node * pre = m_head.load()){
        m_head.store(pre->next);
        delete pre;
    }
}

template <typename Task>
inline std::unique_ptr<Task> LockFreeQueue<Task>::Pop()
{
    CountedNodePtr preHead = m_head.load(std::memory_order_relaxed);
    for(;;){
        IncreaseExternalCount_(m_head, preHead);
        const Node * node = preHead.node;
        if(node == m_tail.load().node){
            node->ReleaseRef();
            return std::unique_ptr<Task>();
        }
        if(m_head.compare_exchange_strong(preHead, node->next)){
            const Task * task = node->task.exchange(nullptr);
            FreeExternalCount_(preHead);
            return std::unique_ptr<Task>(task);
        }
        node->ReleaseRef();
    }
}

template <typename Task>
inline void LockFreeQueue<Task>::Push(Task task)
{
    std::unique_ptr<Task> newTask(new Task(task));
    CountedNodePtr newNext;
    newNext.node = new Node;
    newNext.externalCount = 1;
    CountedNodePtr preTail = m_tail.load();
    for(;;){
        IncreaseExternalCount_(m_tail, preTail);
        Task * preTask = nullptr;
        if(preTail.node->task.compare_exchange_strong(preTask, newTask.get())){
            CountedNodePtr preNext = { 0 };
            if(!preTail.Node->next.compare_exchange_strong(preNext, newNext)){
                delete preNext.node;
                newNext = preNext;
            }
            SetNewTail_(preTail, newNext);
            newTask.release();
            break;
        }
        else{
            CountedNodePtr preNext = { 0 };
            if(preTail.node->next.compare_exchange_strong(preNext, newNext)){
                preNext = newNext;
                newNext.node = new Node;
            }
            SetNewTail_(preTail, preNext);
        }
    }
}

template <typename Task>
inline void LockFreeQueue<Task>::SetNewTail_(CountedNodePtr & preTail, const CountedNodePtr & newTail)
{
    const Node * currTailNode = preTail.node;
    while(!m_tail.compare_exchange_weak(preTail, newTail) && preTail.node == currTailNode);
    if(preTail.node == currTailNode) FreeExternalCount_(preTail);
    else currTailNode->ReleaseRef();
}


template <typename Task>
inline void LockFreeQueue<Task>::IncreaseExternalCount_(std::atomic<CountedNodePtr> & counter, CountedNodePtr & preCounter)
{
    CountedNodePtr newCounter;
    do {
        newCounter = preCounter;
        ++newCounter.externalCount;
    } while(!counter.compare_exchange_strong(preCounter, newCounter, std::memory_order_acquire, std::memory_order_relaxed));
    preCounter.externalCount = newCounter.externalCount;
}

template <typename Task>
inline void LockFreeQueue<Task>::FreeExternalCount_(CountedNodePtr & preNodePtr)
{
    const Node * node = preNodePtr.node;
    const int countIncrease = preNodePtr.externalCount -2;
    NodeCounter preCounter = node->count.load(std::memory_order_relaxed);
    NodeCounter newCounter;
    do {
        newCounter = preCounter;
        --newCounter.externalCounters;
    } while(!node->count.compare_exchange_strong(preCounter, newCounter, std::memory_order_acquire, std::memory_order_relaxed));
    if(!newCounter.internalCount && !newCounter.externalCounters) { delete node; }
}

}//namespace thread
}//namespace generic