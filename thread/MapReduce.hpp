/**
 * @file MapReduce.hpp
 * @author bwu
 * @brief A header only map reduce library on single-machine platform 
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Exception.hpp"
#include "generic/thread/ThreadPool.hpp"
#include <boost/iterator/iterator_facade.hpp>
#include <boost/core/noncopyable.hpp>
#include <boost/functional/hash.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <deque>
#include <map>
namespace generic {
namespace thread {
///@brief namespace of map reduce library
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
    std::string outputFileSpec = "mapreduce_";
    std::string inputFilesDir = "";
    std::streamsize maxFileSegments = 1048576L;

    Specification() {}
};

struct Results
{
    struct Counters
    {
        size_t mapKeysExecuted{0};
        size_t mapKeyErrors{0};
        size_t mapKeysCompleted{0};

        size_t reduceKeysExecuted{0};
        size_t reduceKeyErrors{0};
        size_t reduceKeysCompleted{0};

        size_t numResultFiles = 0;
        Counters(){}
    };

    using Duration = std::chrono::duration<double>;

    Counters counters;
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

struct NullLocker
{
    void lock(){}
    void unlock(){}
    bool try_lock() const { return true; }
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

template <typename Key, typename Value>
class FileHandler
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
struct FileHandler<std::string, std::ifstream>::Data {};

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

template <typename MapTask,
          typename Handle = FileHandler<typename MapTask::Key, typename MapTask::Value> >
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
          typename Partitioner = HashPartitioner,
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
                m_current.first < m_outer->m_partitions && m_iterators[m_current.first] == m_outer->m_intermediates[m_current.first].end();
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
        std::swap(m_intermediates, other.m_intermediates);
    }

    void RunIntermediateResultsShuffle(size_t)
    {
    }

    template <typename Callback>
    void Reduce(size_t partition, Callback & callback)
    {
        typename Intermediates::value_type map;
        std::swap(map, m_intermediates[partition]);

        for(const auto & result : map)
            callback(result.first, result.second.cbegin(), result.second.cend());
    }

    void MergeFrom(size_t partition, InMemory & other)
    {
        auto & thisMap = m_intermediates[partition];
        auto & otherMap = other.m_intermediates[partition];
        if(thisMap.size() == 0){
            std::swap(thisMap, otherMap);
            return;
        }

        using MappedType = typename Intermediates::value_type::mapped_type;
        for(const auto & result : otherMap){
            auto iter = thisMap.insert(std::make_pair(result.first, MappedType())).first;
            std::copy(result.second.cbegin(), result.second.cend(), std::back_inserter(iter->second));
        }
    }

    void MergeFrom(InMemory & other)
    {
        for(size_t i = 0; i < m_partitions; ++i)
            MergeFrom(i, other);
        other.m_intermediates.clear();
    }

    template <typename T>
    bool Insert(const Key & key, const Value & value, T & storeResult)
    {
        return storeResult(key, value) && Insert(key, value);
    }

    bool Insert(const Key & key, const Value & value)
    {
        size_t partition = (m_partitions == 1) ? 0 : m_partitioner(key, m_partitions);
        auto & map = m_intermediates[partition];

        using MappedType = typename Intermediates::value_type::mapped_type;
        auto tmp = std::make_pair(key, MappedType());
        map.insert(std::move(tmp)).first->second.push_back(value);
        return true;
    }

    template <typename FnObj>
    void Combine(FnObj & fnObj)
    {
        Intermediates intermediates;
        intermediates.resize(m_partitions);
        std::swap(m_intermediates, intermediates);

        for(const auto & intermediate : intermediates){
            for(const auto & kv : intermediate){
                fnObj.Start(kv.first);
                for(const auto & value : kv.second)
                    fnObj(value);
                fnObj.Finish(kv.first, *this);
            }
        }
    }

    void Combine(NullCombiner &)
    {
    }

private:
    const size_t m_partitions;
    Partitioner m_partitioner;
    Intermediates m_intermediates;
};


template <typename MapTask, typename ReduceTask>
class ReduceFileOut
{
public:
    ReduceFileOut(const std::string & outputFileSpec, size_t partition, size_t partitions)
    {
        m_filename = outputFileSpec + std::to_string(partition + 1) + "_of_" + std::to_string(partitions);
        m_outputFile.open(m_filename.c_str(), std::ios_base::binary);
        if(!m_outputFile.is_open())
            GENERIC_THROW(std::runtime_error("Error: fail to open file " + m_filename))
    }

    void operator() (typename ReduceTask::Key & key, typename ReduceTask::Value & value)
    {
        m_outputFile << key << "\t" << value << "\r";
    }

private:
    std::string m_filename;
    std::ofstream m_outputFile;
};

template <typename T>
struct SharedPtrLess
{
    bool operator() (const std::shared_ptr<T> & left, const std::shared_ptr<T> & right) const
    {
        return *left < *right;
    }
};

template <typename Record>
inline bool FileKeyCombiner(const std::string & in, const std::string & out, const size_t maxLines = 4294967000U)
{
    std::deque<std::string> tmpFiles;
    
    std::ifstream inFile(in.c_str(), std::ios_base::in | std::ios_base::binary);
    if(!inFile.is_open()){
        GENERIC_THROW(std::runtime_error("Error: fail to open file " + in))
    }
    
    using LinesType = std::map<std::shared_ptr<Record>, std::streamsize, SharedPtrLess<Record> >;
    while(!inFile.eof()){
        LinesType lines;
        for(size_t i = 0; !inFile.eof() && i < maxLines; ++i){
            if(inFile.fail()){
                GENERIC_THROW(std::runtime_error("Error: fail to read file " + in))
            }

            std::string line;
            std::getline(inFile, line);
            if(line.length() > 0){
                //todo
            }
        }

    }

}

// template <typename MapTask, typename ReduceTask,
//           typename KeyType = typename ReduceTask::Key,
//           typename Partitioner = HashPartitioner,
//           typename StoreResultType = ReduceFileOut<MapTask, ReduceTask>,
//           typename 


}//namespace intermediate

/**
 * @brief 
 * 
 * @tparam MapTaskType implements a mapping function to process key/value pairs to generate a set of intermediate key/value pairs
 * @tparam ReduceTaskType implements a reduce function to merge all intermediate values associated with the same intermediate key
 * @note ReduceTaskType should define the key/value types for the results of the reduce phase
 * @tparam DataSourceType implements a mechanism to feed data to the map tasks
 * @tparam CombinerType can by used to partially consolidate results of the map task before they are passed to the reduce tasks
 * @tparam IntermediateStoreType handles storage, merging and sorting of intermediate results between the map and reduce phases
 * @tparam StoreResult intermediate data result type
 */
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
    using StoreResult = StoreResultType;
    using KeyValuePair = typename IntermediateStore::KeyValuePair;
    using ConstResultIterator = typename IntermediateStore::ConstResultIterator;

private:
    class MapTaskRunner : private boost::noncopyable
    {
    public:
        MapTaskRunner(Job & job)
         : m_job(job)
         , m_store(job.NumOfPartitions())
        {
        }

        MapTaskRunner & operator() (const typename MapTask::Key & key, typename MapTask::Value & value)
        {
            MapTask()(*this, key, value);

            Combiner combiner;
            m_store.Combine(combiner);

            return *this;
        }

        template <typename T>
        bool EmitIntermediate(const T & key, const typename ReduceTask::Value & value)
        {
            return m_store.Insert(key, value);
        }

        IntermediateStore & GetIntermediateStore()
        {
            return m_store;
        }

    private:
        Job & m_job;
        IntermediateStore m_store;
    };

    class ReduceTaskRunner : private boost::noncopyable
    {
    public:
        ReduceTaskRunner(std::string outputFileSpec,
                         size_t partition, size_t partitions,
                         IntermediateStore & store, Results & results)
         : m_partition(partition)
         , m_results(results)
         , m_store(store)
         , m_storeResult(outputFileSpec, partition, partitions)
        {
        }

        void Reduce()
        {
            m_store.Reduce(m_partition, *this);
        }

        void Emit(const typename ReduceTask::Key & key, const typename ReduceTask::Value & value)
        {
            m_store.Insert(key, value, m_storeResult);
        }

        template <typename Iterator>
        void operator() (const typename ReduceTask::Key & key, Iterator begin, Iterator end)
        {
            m_results.counters.reduceKeysExecuted++;
            ReduceTask()(*this, key, begin, end);
            m_results.counters.reduceKeysCompleted++;
        }

    private:
        size_t m_partition;
        Results & m_results;
        IntermediateStore & m_store;
        StoreResult m_storeResult;
    };

public:
    Job(DataSource & source, const Specification & spec)
     : m_source(source)
     , m_spec(spec)
     , m_store(IntermediateStore(spec.reduceTasks))
    {
    }

    ConstResultIterator BeginResults() const
    {
        return m_store.BeginResults();
    }

    ConstResultIterator EndResults() const
    {
        return m_store.EndResults();
    }

    bool GetNextMapKey(typename MapTask::Key * & key)
    {
        auto nextKey = std::make_unique<typename MapTask::Key>();
        if(!m_source.SetupKey(*nextKey)) return false;
        key = nextKey.release();
        return true;
    }

    size_t NumOfPartitions() const
    {
        return m_spec.reduceTasks;
    }

    size_t NumOfMapTasks() const
    {
        return m_spec.mapTasks;
    }

    template <typename Schedule>
    void Run(Results & results)
    {
        Schedule schedule;
        Run(schedule, results);
    }

    template <typename Schedule>
    void Run(Schedule && schedule, Results & results)
    {
        const auto startT = std::chrono::system_clock::now();
        schedule(*this, results);
        results.jobRuntime = std::chrono::system_clock::now() - startT;
    }

    template <typename Sync>
    bool RunMapTask(typename MapTask::Key * key, Results & results, Sync & sync)
    {
        const auto startTime = std::chrono::system_clock::now();

        try {
            results.counters.mapKeysExecuted++;

            auto pKey = std::unique_ptr<typename MapTask::Key>(key);
            typename MapTask::Value value;
            if(!m_source.GetData(*pKey, value)){
                results.counters.mapKeyErrors++;
                return false;
            }

            MapTaskRunner runner(*this);
            runner(*pKey, value);

            std::lock_guard<Sync> lock(sync);
            m_store.MergeFrom(runner.GetIntermediateStore());
            results.counters.mapKeysCompleted++;
        }
        catch (std::exception & e){
            std::cerr << "Error: " << e.what() << std::endl;
            results.counters.mapKeyErrors++;
            return false;
        }
        results.mapTimes.push_back(std::chrono::system_clock::now() - startTime);
        return true;
    }

    void RunIntermediateResultsShuffle(size_t partition)
    {
        m_store.RunIntermediateResultsShuffle(partition);
    }

    bool RunReduceTask(size_t partition, Results & results)
    {
        bool success = true;
        const auto & startTime(std::chrono::system_clock::now());
        try {
            ReduceTaskRunner runner(m_spec.outputFileSpec, partition, NumOfPartitions(), m_store, results);
            runner.Reduce();
        }
        catch (std::exception & e){
            std::cerr << "Error: " << e.what() << std::endl;
            results.counters.reduceKeyErrors++;
            success = false;
        }
        results.reduceTimes.push_back(std::chrono::system_clock::now() - startTime);
        return success;
    }

private:
    DataSource & m_source;
    const Specification & m_spec;
    IntermediateStore m_store;
};

namespace schedule {

///@brief Schedule policy that run one map task followed by one reduce task, useful for debug
template <typename Job>
class Sequential
{
public:
    void operator()(Job & job, Results & results)
    {
        Map(job, results);
        Intermediate(job, results);
        Reduce(job, results);
    }

private:
    void Map(Job & job, Results & results)
    {
        const auto startTime(std::chrono::system_clock::now());

        typename Job::MapTask::Key * key = nullptr;
        NullLocker locker;
        while(job.GetNextMapKey(key) && job.RunMapTask(key, results, locker)){}
        results.mapRuntime = std::chrono::system_clock::now() - startTime;
    }

    void Intermediate(Job & job, Results & results)
    {
        const auto startTime(std::chrono::system_clock::now());
        for(size_t i = 0; i < job.NumOfPartitions(); ++i)
            job.RunIntermediateResultsShuffle(i);
        results.shuffleRuntime = std::chrono::system_clock::now() - startTime;
    }

    void Reduce(Job & job, Results & results)
    {
        const auto startTime(std::chrono::system_clock::now());
        for(size_t i = 0; i < job.NumOfPartitions(); ++i)
            job.RunReduceTask(i, results);
        results.reduceRuntime = std::chrono::system_clock::now() - startTime;
    }
};

///@brief Schedule policy that uses avaliable cores to run as many map simultaneous tasks as possible
template <typename Job>
class Parallel
{
    using LocalResults = std::vector<std::shared_ptr<Results> >;
    size_t m_threads;
    LocalResults m_localResults;
public:
    Parallel()
    {
        m_threads = std::thread::hardware_concurrency();
        if(m_threads > 1) m_threads -= 1;
    }

    void SetThreads(size_t threads)
    {
        m_threads = std::min<size_t>(std::thread::hardware_concurrency(), threads);
        if(0 == m_threads) m_threads = 1;
    }

    void operator()(Job & job, Results & results)
    {
        Map(job, results);
        Intermediate(job, results);
        Reduce(job, results);
        CollateResults(results);
        results.counters.numResultFiles = job.NumOfPartitions();
    }

private:
    void Map(Job & job, Results & results)
    {
        auto startTime = std::chrono::system_clock::now();
        size_t threads = std::min(m_threads, job.NumOfMapTasks());

        std::mutex lock;
        ThreadPool pool(threads);
        typename Job::MapTask::Key * key = nullptr;
        while(job.GetNextMapKey(key)){
            pool.Submit(std::bind(&Job::template RunMapTask<std::mutex>, std::ref(job), key, std::ref(results), std::ref(lock)));
        }
        pool.Wait();
        results.mapRuntime = std::chrono::system_clock::now() - startTime;
    }

    void Intermediate(Job & job, Results & results)
    {
        const auto startTime(std::chrono::system_clock::now());
        size_t threads = std::min(m_threads, job.NumOfPartitions());

        ThreadPool pool(threads);
        for(size_t i = 0; i < job.NumOfPartitions(); ++i){
            pool.Submit(std::bind(&Job::RunIntermediateResultsShuffle, std::ref(job), i));
        }
        pool.Wait();
        results.shuffleRuntime = std::chrono::system_clock::now() - startTime;
    }

    void Reduce(Job & job, Results & results)
    {
        const auto startTime(std::chrono::system_clock::now());
        size_t threads = std::min(m_threads, job.NumOfPartitions());

        ThreadPool pool(threads);
        for(size_t i = 0; i < job.NumOfPartitions(); ++i){
            auto localResult = std::make_shared<Results>();
            m_localResults.push_back(localResult);
            pool.Submit(std::bind(&Job::RunReduceTask, std::ref(job), i, std::ref(*localResult)));
        }
        pool.Wait();
        results.reduceRuntime = std::chrono::system_clock::now() - startTime;
    }

    void CollateResults(Results & results)
    {
        for(const auto & localRes : m_localResults){
            results.counters.mapKeysCompleted += localRes->counters.mapKeysCompleted;
            results.counters.mapKeysExecuted += localRes->counters.mapKeysExecuted;
            results.counters.mapKeyErrors += localRes->counters.mapKeyErrors;

            results.counters.reduceKeysCompleted += localRes->counters.reduceKeysCompleted;
            results.counters.reduceKeysExecuted += localRes->counters.reduceKeysExecuted;
            results.counters.reduceKeyErrors += localRes->counters.reduceKeyErrors;

            std::copy(localRes->mapTimes.cbegin(), localRes->mapTimes.cend(), std::back_inserter(results.mapTimes));
            std::copy(localRes->shuffleTimes.cbegin(), localRes->shuffleTimes.cend(), std::back_inserter(results.shuffleTimes));
            std::copy(localRes->reduceTimes.cbegin(), localRes->reduceTimes.cend(), std::back_inserter(results.reduceTimes));
        }
        m_localResults.clear();
    }
};

}//namespace schedule

}//mapreduce
}//thread
}//generic