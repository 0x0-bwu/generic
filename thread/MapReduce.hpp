#ifndef GENERIC_THREAD_MAPREDUCE_MAPREDUCE_HPP
#define GENERIC_THREAD_MAPREDUCE_MAPREDUCE_HPP
#include "generic/common/Exception.hpp"
#include "generic/thread/ThreadPool.hpp"
#include <boost/iterator/iterator_facade.hpp>
#include <boost/core/noncopyable.hpp>
#include <boost/functional/hash.hpp>
#include <fstream>
#include <chrono>
#include <vector>
#include <map>
namespace generic {
namespace thread {
namespace mapreduce {

struct HashPartitioner
{
    template <typename T>
    size_t operator()(const T & key, size_t partitions) const
    {
        boost::hash<T> hasher;
        return hasher(key) % partitions;
    }
};
struct Specification
{
    size_t mapTasks = 0;    // ideal number of map tasks to use
    size_t reduceTasks = 1; // ideal number of reduce tasks to use

    Specification() {}
};

struct Results
{
    struct TagCounters
    {
        size_t actualMapTasks = 0;
        size_t actualReduceTasks = 0;

        size_t mapKeysExecuted = 0;
        size_t mapKeyErrors = 0;
        size_t mapKeysCompleted = 0;

        size_t reduceKeysExecuted = 0;
        size_t reduceKeyErrors = 0;
        size_t reduceKeysCompleted = 0;

        TagCounters(){}
    };

    using Duration = std::chrono::duration<double>;

    TagCounters counters;
    Duration jobRuntime;
    Duration mapRuntime;
    Duration shuffleRuntime;
    Duration reduceRuntime;
    std::vector<Duration> mapTimes;
    std::vector<Duration> shuffleTimes;
    std::vector<Duration> reduceTimes;
};

template <typename KeyType, typename ValueType>
struct MapTask
{
    using Key = KeyType;
    using Value = ValueType;
};

template <typename KeyType, typename ValueType>
struct ReduceTask
{
    using Key = KeyType;
    using Value = ValueType;
};

struct NullCombiner
{
    template <typename IntermediateStore>
    static void Run(IntermediateStore &){}

    template <typename ReduceTaskKey>
    void Start(const ReduceTaskKey &){}

    template <typename ReduceTaskKey, typename IntermediateStore>
    void Finish(const ReduceTaskKey &, IntermediateStore &){}

    template <typename ReduceTaskKey>
    void operator() (const ReduceTaskKey &){}
};

namespace datasource {
namespace detail {
template <typename Key, typename Value>
class FileHandler : private boost::nocopyable
{
public:
    FileHandler(const Specification & spec);

    bool GetData(const Key & key, Value & value) const;
    bool SetupKey(Key &) const { return false; }

private:
    const Specification & m_spec;

    struct Data;
    std::shared_ptr<Data> m_data;
};

template <>
inline struct FileHandler<std::string, std::ifstream>::Data {};

template <>
inline FileHandler<std::string, std::ifstream>::FileHandler(const Specification & spec)
 : m_spec(spec)
 , m_data(new Data)
{
}

template <>
inline bool FileHandler<std::string, std::ifstream>::GetData(const std::string & key, std::ifstream & value) const
{
    value.open(key.c_str(), std::ios_base::binary);
    return value.is_open();
}
}//namespace detail

template <typename MapTask,
          typename Handle = detail::FileHandle<typename MapTask::Key, typename MapTask::Value> >
class DirectoryIterator : private boost::noncopyable
{
public:
    DirectoryIterator(Specification & spec){}//todo
};

}//namespace datasource

namespace intermediate
{

template <typename ReduceKeyType, typename MapValueType>
inline ReduceKeyType makeIntermediateKey(const MapValueType &);

template <>
inline std::string makeIntermediateKey(const std::pair<const char *, std::uintmax_t> & value)
{
    GENERIC_ASSERT(value.second < std::numeric_limits<std::string::size_type>::max())
    return std::string(value.first, (std::string::size_type)value.second);
}

template <>
inline std::pair<char const *, std::uintmax_t> makeIntermediateKey(const std::string & value)
{
    return std::make_pair(value.c_str(), value.length());
}

template <typename MapTask, typename ReduceTask>
class ReduceNullOutput
{
public:
    ReduceNullOutput(const std::string &, size_t, size_t){}

    bool operator() (const typename MapTask::Key &, const typename ReduceTask::Value &)
    {
        return true;
    }
};

template <typename MapTask, typename ReduceTask,
          typename KeyType = typename ReduceTask::Key,
          typename PartitionFunc = HashPartitioner,
          typename KeyCompare = std::less<typename ReduceTask::Key>,
          typename StoreResultType = ReduceNullOutput<MapTask, ReduceTask> >
class InMemory : private boost::noncopyable
{
public:
    using Key = KeyType;
    using Value = typename ReduceTask::Value;
    using KeyValuePair = std::pair<Key, Value>;
    using StoreResult = StoreResultType;
private:
    using Intermediate = std::map<Key, std::list<Value>, KeyCompare>;
    using Intermediates = std::vector<Intermediate>;

public:
    class ConstResultIterator
     : public boost::iterator_facade<ConstResultIterator, const KeyValuePair, boost::forward_traversal_tag>
    {
        friend class InMemory;
        friend class boost::iterator_core_access;
    public:
        ConstResultIterator(const ConstResultIterator &) = default;
        ConstResultIterator & operator= (const ConstResultIterator &) = default;
    
    private:
        explicit ConstResultIterator(const InMemory * outer)
         : m_outer(outer)
        {
            GENERIC_ASSERT(outer)
            m_iterators.resize(outer->m_partitions);
        }

        void increment()
        {
            ++m_current.second;
            if(m_current.second == m_iterators[m_current.first]->second.end()){
                if(m_iterators[m_current.first] != m_outer->m_intermediates[m_current.first].end())
                    ++m_iterators[m_current.first];
                set_current();
            }
            else
                m_value = std::make_pair(m_iterators[m_current.first]->first, *m_current.second);
        }

        bool equal(const ConstResultIterator & other) const
        {
            if(m_current.first == std::numeric_limits<size_t>::max() ||
                other.m_current.first == std::numeric_limits<size_t>::max())
                return other.m_current.first == m_current.first;
            return m_value == other.m_value;
        }

        ConstResultIterator & begin()
        {
            for(size_t loop = 0; loop < m_outer->m_partitions; ++loop)
                m_iterators[loop] = m_outer->m_intermediates[loop].cbegin();
            set_current();
            return *this;
        }

        ConstResultIterator & end()
        {
            m_current.first = std::numeric_limits<size_t>::max();
            m_value = KeyValuePair{};
            m_iterators.clear();
            return *this;
        }

        const KeyValuePair & dereference() const
        {
            return m_value;
        }

        void set_current()
        {
            for(m_current.first = 0;
                m_current.first < m_outer->m_partitions && m_iterators[m_current.first] == m_outer->m_intermediates[m_current.first].end());
                ++m_current.first){}
            
            for(auto loop = m_current.first + 1; loop < m_outer->m_partitions; ++loop){
                if(m_iterators[loop] != m_outer->m_intermediates[loop].end() &&
                    *m_iterators[m_current.first] > *m_iterators[loop])
                    m_current.first = loop;
            }

            if(m_current.first == m_outer->m_partitions) end();
            else{
                m_current.second = m_iterators[m_current.first]->second.cbegin();
                m_value = std::make_pair(m_iterators[m_current.first]->first, *m_current.second);
            }
        }

    private:
        using Iterators = std::vector<typename Intermediates::value_type::const_iterator>;
        using Current = std::pair<size_t, typename Intermediates::value_type::mapped_type::const_iterator>;

        KeyValuePair m_value;
        Iterators m_iterators;
        const InMemory * m_outer;
        Current m_current;
    };
    friend class ConstResultIterator;

    explicit InMemory(size_t partitions = 1)
     : m_partitions(partitions)
    {
        m_intermediates.resize(m_partitions);
    }

    ConstResultIterator BeginResults() const
    {
        return ConstResultIterator(this).begin();
    }

    ConstResultIterator EndResults() const
    {
        return ConstResultIterator(this).end();
    }

    void Swap(InMemory & other)
    {
        swap(m_intermediates, other.m_intermediates);
    }

    void RunIntermediateResultsShuffle(size_t)
    {
    }
    //todo


private:
    const size_t m_partitions;
    Intermediates m_intermediates;
};

}//namespace intermediate

template <typename MapTaskType, 
          typename ReduceTaskType,
          typename DataSourceType,
          typename CombinerType = NullCombiner,
          typename IntermediateStoreType = intermediate::InMemory<MapTaskType, ReduceTaskType>,
          typename StoreResultType = typename IntermediateStoreType::StoreResult>
class Job : private boost::noncopyable
{
public:
    using MapTask = MapTaskType;
    using ReduceTask = ReduceTaskType;
    using DataSource = DataSourceType;
    using Combiner = CombinerType;
    using IntermediateStore = IntermediateStoreType;
private:
    class MapTaskRunner : private boost::noncopyable
    {
    public:

    };

public:
    
}


namespace schedule {

template <typename Job>
class Sequential
{
public:
    void operator()(Job & job, Results & results)
    {
        //todo
    }

    void Map(Job & job, Results & results)
    {
        const auto startTime(std::chrono::system_clock::now());

        //todo
    }
};

}//namespace schedule

}//mapreduce
}//thread
}//generic

#endif//GENERIC_THREAD_MAPREDUCE_MAPREDUCE_HPP