/**
 * @file ZipView.hpp
 * @author bwu
 * @brief simple implementation of c++23's std::ranges::zip
 * @version 0.1
 * @date 2025-01-08
 */
#pragma once
#include <tuple>
#include <vector>
#include <iterator>
namespace generic::utils {
namespace impl {

template <typename... Iterators>
class ZipIterator
{
public:
    using ValueType = std::tuple<typename std::iterator_traits<Iterators>::reference...>;
/**
 * @brief Brief description of ZipIterator.
 * @param m_iterators(iterators...
 * @return explicit
 */
    explicit ZipIterator(Iterators... iterators) : m_iterators(iterators...) {}

    ZipIterator & operator ++ ()
    {
        Increment(std::index_sequence_for<Iterators...>{});
        return *this;
    }

    bool operator != (const ZipIterator & other) const
    {
        return NE(other, std::index_sequence_for<Iterators...>{});
    }

    ValueType operator * () const
    {
        return Dereference(std::index_sequence_for<Iterators...>{});
    }

private:
    template <std::size_t... Is>
    void Increment(std::index_sequence<Is...>)
    {
        (++std::get<Is>(m_iterators), ...);
    }

    template <std::size_t... Is>
    bool NE(const ZipIterator & other, std::index_sequence<Is...>) const
    {
        return (... || (std::get<Is>(m_iterators) != std::get<Is>(other.m_iterators)));
    }

    template <std::size_t... Is>
    ValueType Dereference(std::index_sequence<Is...>) const {
        return std::tie(*std::get<Is>(m_iterators)...);
    }
private:
    std::tuple<Iterators...> m_iterators;
};

template <typename... Containers>
class ZipView {
public:
    using Iterator = ZipIterator<decltype(std::declval<Containers&>().begin())...>;

    ZipView(Containers &... containers) : m_containers(containers...) {}

    Iterator begin()
    {
        return BeginImpl(std::index_sequence_for<Containers...>{});
    }

    Iterator end()
    {
        return EndImpl(std::index_sequence_for<Containers...>{});
    }

private:
    template <std::size_t... Is>
    Iterator BeginImpl(std::index_sequence<Is...>)
    {
        return Iterator(std::get<Is>(m_containers).begin() ...);
    }

    template <std::size_t... Is>
    Iterator EndImpl(std::index_sequence<Is...>)
    {
        return Iterator(std::get<Is>(m_containers).end() ...);
    }

    std::tuple<Containers&...> m_containers;
};

} // namespace impl

template <typename... Containers>
impl::ZipView<Containers...> Zip(Containers & ... containers)
{
    return impl::ZipView<Containers...>(containers...);
}

} // namespace generic::utils

