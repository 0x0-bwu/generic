//for test
#include "thread/ThreadPool.hpp"
#include "math/MathUtility.hpp"
#include <unordered_map>
#include <shared_mutex>
#include <iostream>
using namespace generic;

class Cache
{
public:
    std::string GetOrCreate(size_t i) {
        {
            std::shared_lock<std::shared_mutex> lock(m_mutex);
            auto iter = m_cache.find(i);
            if (iter != m_cache.end())
                return iter->second;
        }
        return create(i);
    }

    std::string create(size_t i) {
        std::unique_lock<std::shared_mutex> lock(m_mutex);
        auto res = m_cache.emplace(i, std::to_string(i));
        return res.first->second;
    }

    size_t size() const {
        // std::shared_lock<std::shared_mutex> lock(m_mutex);
        return m_cache.size();
    }

private:
    std::shared_mutex m_mutex;
    std::unordered_map<size_t, std::string> m_cache;

};

int main(int argc, char * argv[])
{
    //test code here
    thread::ThreadPool pool(4);

    Cache cache;
    auto task = [&cache]{
        for(size_t i = 0; i < 1000000; ++i) {
            cache.GetOrCreate(math::Random(0, 20000));
        }
    };

    for(size_t i = 0; i < 6; ++i) {
        pool.Submit(task);
    }

    pool.Wait();
    std::cout << cache.size() << std::endl;
}