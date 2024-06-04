/**
 * @file RandomTree.hpp
 * @author bwu
 * @brief Model of random tree
 * @version 0.1
 * @date 2024-05-20
 */
#pragma once
#include "generic/math/MathUtility.hpp"
#include "generic/common/Traits.hpp"
#include <unordered_set>
#include <queue>

namespace generic::tree {

class RandomTree
{
public:
    struct TagSilkworm{};
    struct TagRandom{};
    struct TagChain{};
    struct TagTall{};
    std::vector<size_t> p, nid;
    std::unordered_set<size_t> leaves;
    inline static constexpr size_t invalid = std::numeric_limits<size_t>::max();
    RandomTree()
    {
        Init(1);
    }

    explicit RandomTree(size_t n, [[maybe_unused]] TagRandom t = TagRandom{})
    {
        Init(n); 
        Random(n - 1, 0);     
    }

    explicit RandomTree(size_t n, size_t k, TagTall)
    {
        Init(n);
        Tall(n - 1, k, 0);
    }

    explicit RandomTree(size_t n, TagChain)
    {
        Init(n); 
        Chain(n - 1, 0);
    }

    explicit RandomTree(size_t n, TagSilkworm)
    {
        Init(n);
        Silkworm(n - 1, 0);
    }

    void Init(size_t n)
    {
        p.reserve(n);
        nid.reserve(n);
        nid.emplace_back(0);
        leaves.emplace(0);
        p.emplace_back(invalid);
    }

    size_t Size() const { return nid.size(); }

    void AddNode(size_t pa)
    {
        GENERIC_ASSERT(pa < Size())
        if (auto iter = leaves.find(pa); iter != leaves.end())
            leaves.erase(iter);
        leaves.emplace(nid.size());
        nid.emplace_back(nid.size());
        p.emplace_back(pa);
    }

    void Random(size_t n, size_t pa)
    {
        size_t sz = Size();
        GENERIC_ASSERT(pa < sz)
        AddNode(pa);
        if (n == 1) return;
        if (n == 2) {
            AddNode(sz);
            return;
        }

        std::vector<size_t> prufer, cnt;
        std::vector<std::vector<size_t> > g;
        g.resize(n);
        cnt.assign(n, 0);
        for (size_t i = 0; i < n - 2; ++i) {
            auto x = math::Random<size_t>(0, n - 1);
            prufer.emplace_back(x);
            ++cnt[x];
        }
        std::priority_queue<size_t> q;
        for (size_t i = 0; i < n; ++i)
            if (not cnt.at(i)) q.push(i);
        for (auto v : prufer) {
            auto u = q.top();
            g[u].emplace_back(v);
            g[v].emplace_back(u);
            q.pop();
            if (--cnt[v] == 0) q.push(v);
        }
        size_t x = q.top();
        q.pop();
        size_t y = q.top();
        g[x].emplace_back(y);
        g[y].emplace_back(x);

        std::queue<size_t> bfs; bfs.push(0);
        size_t id = sz;
        
        while (not bfs.empty()) {
            auto u = bfs.front();
            cnt[u] = 1;
            bfs.pop();
            for (auto v : g[u]) {
                if (cnt[v] == 0) {
                    AddNode(id);
                    bfs.push(v);
                }
            }
            ++id;
        }
    }

    void Tall(size_t n, size_t k, size_t pa)
    {
        auto sz = Size();
        AddNode(pa);
        for (size_t i = sz + 1; i < sz + n; ++i)
            AddNode(math::Random(std::max(sz, i - k), i - 1));
    }

    void Chain(size_t n, size_t pa)
    {
        Tall(n, 1, pa);
    }

    void Silkworm(size_t n, size_t pa)
    {
        auto sz = Size();
        auto len = (n + 1) / 2;
        Chain(len, pa);
        for (size_t i = sz; i + len < sz + n; ++i)
            AddNode(i);
    }
};

} // generic::tree