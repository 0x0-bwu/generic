/**
 * @file LockFreeBitSet.hpp
 * @author bwu
 * @brief Lock free bit set, modified from folly concurrent bit set
 * @version 0.1
 * @date 2024-11-19
 */
#pragma once
#include "generic/common/Exception.hpp"
#include <limits>
#include <atomic>
#include <array>
namespace generic::thread {

template <size_t N, typename Block = uint64_t>
class LockFreeBitSet
{
	static_assert(std::atomic<Block>::is_always_lock_free/*shoule be lock free type*/);
	static_assert(std::is_integral_v<Block> and std::is_unsigned_v<Block>/*shoule be unsigned integral type*/);
public:
	LockFreeBitSet() { std::fill(m_data.begin(), m_data.end(), 0); }
	LockFreeBitSet(const LockFreeBitSet &) = delete;
	LockFreeBitSet & operator=(const LockFreeBitSet &) = delete;

	/**
	 * Set bit idx to true, using the given memory order. Returns the previous value of the bit.
	 * Note that the operation is a read-modify-write operation due to the use of fetch_or.
	 */
	bool set(size_t idx, std::memory_order order = std::memory_order_seq_cst)
	{
		GENERIC_ASSERT(idx < N);
		Block mask = ONE << BitOffset(idx);
		return m_data[BlockIndex(idx)].fetch_or(mask, order) & mask;
	}

	/**
	 * Set bit idx to false, using the given memory order. Returns the previous value of the bit.
	 * Note that the operation is a read-modify-write operation due to the use of fetch_and.
	 */
	bool reset(size_t idx, std::memory_order order = std::memory_order_seq_cst)
	{
		GENERIC_ASSERT(idx < N);
		Block mask = ONE << BitOffset(idx);
		return m_data[BlockIndex(idx)].fetch_and(~mask, order) & mask;
	}

	/**
	 * Set bit idx to the given value, using the given memory order. Returns the previous value of the bit.
	 * Note that the operation is a read-modify-write operation due to the use of fetch_and or fetch_or.
	 * Yes, this is an overload of set(), to keep as close to std::bitset's interface as possible.
	 */
	bool set(size_t idx, bool value, std::memory_order order = std::memory_order_seq_cst)
	{
		return value ? set(idx, order) : reset(idx, order);
	}

	/// Read bit idx.
	bool test(size_t idx, std::memory_order order = std::memory_order_seq_cst) const
	{
		GENERIC_ASSERT(idx < N);
		Block mask = ONE << BitOffset(idx);
		return m_data[BlockIndex(idx)].load(order) & mask;
	}

	/// Same as test() with the default memory order.
	bool operator[] (size_t idx) const { return test(idx); }

	/// Return the size of the bitset.
	constexpr size_t size() const { return N; }

private:
	static constexpr size_t BITS_PER_BLOCK = std::numeric_limits<Block>::digits;
/**
 * @brief Brief description of BlockIndex.
 * @param bit
 * @return static constexpr size_t
 */
	static constexpr size_t BlockIndex(size_t bit) { return bit / BITS_PER_BLOCK; }
/**
 * @brief Brief description of BitOffset.
 * @param bit
 * @return static constexpr size_t
 */
	static constexpr size_t BitOffset(size_t bit) { return bit % BITS_PER_BLOCK; }

	// avoid casts
	static constexpr Block ONE = 1;
	static constexpr size_t NUM_BLOCKS = (N + BITS_PER_BLOCK - 1) / BITS_PER_BLOCK;
	std::array<std::atomic<Block>, NUM_BLOCKS> m_data;
};

} // namespace generic::thread
