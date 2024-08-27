/**
 * @file Utility.hpp
 * @author bwu
 * @brief Thread related utilities
 * @version 0.1
 * @date 2024-08-26
 */
#pragma once
#include <unordered_map>
#include <shared_mutex>
#include <thread>
#include <memory>
namespace generic::thread {

class FunctionWrapper
{
    struct ImplBase{
        virtual void Call() = 0;
        virtual ~ImplBase(){}
    };

    template <typename Func>
    struct ImplType : ImplBase
    {
        Func fun;
        ImplType(Func && _fun) : fun(std::move(_fun)) {}
        void Call() { fun(); }
    };

public:
    FunctionWrapper() = default;
    FunctionWrapper(FunctionWrapper & other) = delete;
    FunctionWrapper(const FunctionWrapper & other) = delete;
    FunctionWrapper & operator= (const FunctionWrapper & other) = delete;

    FunctionWrapper(FunctionWrapper && other) : m_impl(std::move(other.m_impl)) {}
    FunctionWrapper & operator= (FunctionWrapper && other) { m_impl = std::move(other.m_impl); return *this; }
    
    template <typename F>
    FunctionWrapper(F && f) : m_impl(new ImplType<F>(std::move(f))) {}

    void operator()() { m_impl->Call(); }

private:
    std::unique_ptr<ImplBase> m_impl;
};

class ThreadIdMgr
{
public:
    static ThreadIdMgr & Inst()
    {
        static ThreadIdMgr mgr;
        return mgr;
    }

    static size_t GetId()
    {
        return Inst().GetLocalId();
    }

private:
    size_t GetLocalId()
    {
        {
            std::shared_lock lock(m_mutex);
            if (auto iter = m_localIds.find(std::this_thread::get_id()); iter != m_localIds.cend())
                return iter->second;
        }
        return CreateLocalId();
    }

    size_t CreateLocalId()
    {
        std::unique_lock lock(m_mutex);
        return m_localIds.emplace(std::this_thread::get_id(), m_localIds.size()).first->second;
    }
private:
    ThreadIdMgr() = default;
    std::shared_mutex m_mutex;
    std::unordered_map<std::thread::id, size_t> m_localIds;
};

}//namespace generic::thread
