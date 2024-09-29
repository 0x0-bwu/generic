/**
 * @file LockFree.hpp
 * @author bwu
 * @brief Lock free hash map, modified from folly atomic hash map
 * @version 0.1
 * @date 2024-09-12
 */

#include "generic/common/Exception.hpp"
#include "generic/math/MathUtility.hpp"
#include "Utility.hpp"
#include <memory>
#include <atomic>

namespace generic::thread {

struct LockFreeHashArrayLinearProbeFcn
{
	inline size_t operator() (size_t idx, size_t /* numProbes */, size_t capacity) const
	{
		idx += 1; // linear probing
		// Avoid modulus because it's slow
		return GENERIC_LIKELY(idx < capacity) ? idx : (idx - capacity);
	}
};

struct LockFreeHashArrayQuadraticProbeFcn
{
	inline size_t operator()(size_t idx, size_t numProbes, size_t capacity) const 
	{
		idx += numProbes; // quadratic probing
		// Avoid modulus because it's slow
		return GENERIC_LIKELY(idx < capacity) ? idx : (idx - capacity);
	}
};

namespace detail {

template <typename NotKeyT, typename KeyT>
inline void CheckLegalKeyIfKeyTImpl(NotKeyT /* ignored */, KeyT /* emptyKey */, KeyT /* lockedKey */, KeyT /* erasedKey */) {}

template <typename KeyT>
inline void CheckLegalKeyIfKeyTImpl([[maybe_unused]] KeyT keyIn, [[maybe_unused]] KeyT emptyKey, 
									[[maybe_unused]] KeyT lockedKey, [[maybe_unused]] KeyT erasedKey)
{
  GENERIC_ASSERT(keyIn != emptyKey);
  GENERIC_ASSERT(keyIn != lockedKey);
  GENERIC_ASSERT(keyIn != erasedKey);
}

/**
 * Currently this only supports forward and bidirectional iteration.  The
 * derived class must must have definitions for these methods:
 *
 *   void Increment();
 *   void Decrement(); // optional, to be used with bidirectional
 *   reference Dereference() const;
 *   bool Equal(const [appropriate iterator type] & rhs) const;
 *
 * These names are consistent with those used by the Boost iterator
 * facade / adaptor classes to ease migration classes in this file.
 *
 * Template parameters:
 * D: the deriving class (CRTP)
 * V: value type
 * Tag: the iterator category, one of:
 *   std::forward_iterator_tag
 *   std::bidirectional_iterator_tag
 */
template <class D, class V, class Tag>
class IteratorFacade {
public:
	using value_type = V;
	using reference = value_type&;
	using pointer = value_type*;
	using difference_type = std::ptrdiff_t;
	using iterator_category = Tag;

	friend bool operator==(const D & lhs, const D & rhs) { return Equal(lhs, rhs); }

	friend bool operator!=(const D & lhs, const D & rhs) { return not (lhs == rhs); }

	V & operator* () const { return AsDerivedConst().Dereference(); }

	V * operator-> () const { return std::addressof(operator*()); }

	D & operator++()
	{
		AsDerived().Increment();
		return AsDerived();
	}

	D operator++ (int)
	{
		auto ret = AsDerived(); // copy
		AsDerived().Increment();
		return ret;
	}

	D & operator-- ()
	{
		AsDerived().Decrement();
		return AsDerived();
	}

	D operator--(int)
	{
		auto ret = AsDerived(); // copy
		AsDerived().Decrement();
		return ret;
	}

private:
	D & AsDerived() { return static_cast<D &>(*this); }

	const D & AsDerivedConst() const { return static_cast<const D &>(*this); }

	static bool Equal(const D & lhs, const D & rhs) { return lhs.equal(rhs); }
};

struct Identity
{
  template <class T>
  constexpr T && operator()(T && t) const noexcept { return static_cast<T &&>(t); }
};

} // namespace detail


template <
    class KeyT,
    class ValueT,
    class HashFcn = std::hash<KeyT>,
    class EqualFcn = std::equal_to<KeyT>,
    class Allocator = std::allocator<char>,
    class ProbeFcn = LockFreeHashArrayLinearProbeFcn,
    class KeyConvertFcn = detail::Identity>
class LockFreeHashMap;

/**
 *  LockFreeHashArray is the building block for LockFreeHashMap.  It provides the
 *  core lock-free functionality, but is limited by the fact that it cannot
 *  grow past its initialization size and is a little more awkward (no public
 *  constructor, for example).  If you're confident that you won't run out of
 *  space, don't mind the awkardness, and really need bare-metal performance,
 *  feel free to use LFHA directly.
 */

template <
    class KeyT,
    class ValueT,
    class HashFcn = std::hash<KeyT>,
    class EqualFcn = std::equal_to<KeyT>,
    class Allocator = std::allocator<char>,
    class ProbeFcn = LockFreeHashArrayLinearProbeFcn,
    class KeyConvertFcn = detail::Identity>
class LockFreeHashArray
{
	static_assert((std::is_convertible_v<KeyT, int32_t> or
				   std::is_convertible_v<KeyT, int64_t> or
				   std::is_convertible_v<KeyT, const void*>),
				"You are trying to use LockFreeHashArray with disallowed key types."
				"You must use atomically compare-and-swappable integer keys, or a different container class.");

public:
	using key_type = KeyT;
	using mapped_type = ValueT ;
	using hasher = HashFcn ;
	using key_equal = EqualFcn;
	using key_convert = KeyConvertFcn;
	using value_type = std::pair<const KeyT, ValueT>;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;
	using reference = value_type &;
	using const_reference = const value_type &;
	using pointer = value_type *;
	using const_pointer = const value_type *;

	const size_t capacity;
	const size_t maxEntries;
	const KeyT kEmptyKey;
	const KeyT kLockedKey;
	const KeyT kErasedKey;

	template <class ContT, class IterVal>
	struct Iterator : detail::IteratorFacade<Iterator<ContT, IterVal>, IterVal, std::forward_iterator_tag>
	{
	public:
		Iterator() = default;
		// Conversion ctor for interoperability between const_iterator and
		// iterator.  The enable_if<> magic keeps us well-behaved for
		// is_convertible<> (v. the iterator_facade documentation).
		template <class OtherContT, class OtherVal>
		Iterator(const Iterator<OtherContT, OtherVal> & o, std::enable_if_t<std::is_convertible_v<OtherVal*, IterVal*>> * = nullptr)
	     : m_container(o.m_container), m_offset(o.m_offset)
		{
		}

		explicit Iterator(ContT * array, size_t offset)
	     : m_container(array), m_offset(offset)
		{
		}

		// Returns unique index that can be used with FindAt().
		// WARNING: The following function will fail silently for hashtable
		// with capacity > 2^32
		uint32_t GetIndex() const { return m_offset; }

		void AdvancePastEmpty()
		{
			while (m_offset < m_container->capacity && not isValid()) {
				++m_offset;
			}
		}
	private:
		friend class LockFreeHashArray;
		friend class detail::IteratorFacade<Iterator, IterVal, std::forward_iterator_tag>;

		void Increment()
		{
			++m_offset;
			AdvancePastEmpty();
		}

		bool Equal(const Iterator & other) const
		{
			return m_container == other.m_container && m_offset == other.m_offset;
		}

		IterVal & Dereference() const { return m_container->m_cells[m_offset]; }

		bool isValid() const 
		{
			KeyT key = AcquireLoadKey(m_container->m_cells[m_offset]);
			return key != m_container->kEmptyKey && key != m_container->kLockedKey && key != m_container->kErasedKey;
		}

	private:
		ContT * m_container{nullptr};
		size_t  m_offset{0};
	};

	using const_iterator = Iterator<const LockFreeHashArray, const value_type> ;
	using iterator = Iterator<LockFreeHashArray, value_type>;

	LockFreeHashArray(size_t capacity, KeyT emptyKey, KeyT lockedKey, KeyT erasedKey, double maxLoadFactor, uint32_t cacheSize)
	 : capacity(capacity), maxEntries(size_t(maxLoadFactor * capacity + 0.5)), kEmptyKey(emptyKey),
       kLockedKey(lockedKey), kErasedKey(erasedKey), m_kAnchorMask(math::NextPowTwo(capacity) - 1)
	{
		if (0 == capacity) ThrowException("zero capacity");
	}

	LockFreeHashArray(const LockFreeHashArray &) = delete;
	LockFreeHashArray& operator=(const LockFreeHashArray &) = delete;
	~LockFreeHashArray() = default;

	// You really shouldn't need this if you use the SmartPtr provided by create,
	// but if you really want to do something crazy like stick the released
	// pointer into a DescriminatedPtr or something, you'll need this to clean up
	// after yourself.
  	static void Destroy(LockFreeHashArray * p)
	{
		GENERIC_ASSERT(p);
		size_t sz = sizeof(LockFreeHashArray) + sizeof(value_type) * p->capacity;

		for (size_t i = 0; i < p->capacity; ++i) {
			if (p->m_cells[i].first != p->kEmptyKey) {
				p->m_cells[i].~value_type();
			}
		}
		p->~LockFreeHashArray();

		Allocator().deallocate((char*)p, sz);
	}
private:
 	const size_t m_kAnchorMask;

	struct Deleter
	{
		void operator()(LockFreeHashArray * ptr) { LockFreeHashArray::Destroy(ptr); }
	};

	struct SimpleRetT
	{
		size_t idx{0};
		bool success{false};
		SimpleRetT(size_t i, bool s) : idx(i), success(s) {}
		SimpleRetT() = default;
	};

	static std::atomic<KeyT> * CellKeyPtr(const value_type & r)
	{
		// We need some illegal casting here in order to actually store
		// our value_type as a std::pair<const,>.  But a little bit of
		// undefined behavior never hurt anyone ...
		static_assert(sizeof(std::atomic<KeyT>) == sizeof(KeyT), "std::atomic is implemented in an unexpected way for LFHM");
		return const_cast<std::atomic<KeyT>*>(reinterpret_cast<const std::atomic<KeyT> *>(&r.first));
	}

public:
 	using SmartPtr = std::unique_ptr<LockFreeHashArray, Deleter>;

	/*
	* create --
	*
	*   Creates LockFreeHashArray objects.  Use instead of constructor/destructor.
	*
	*   We do things this way in order to avoid the perf penalty of a second
	*   pointer indirection when composing these into LockFreeHashMap, which needs
	*   to store an array of pointers so that it can perform atomic operations on
	*   them when growing.
	*
	*   Instead of a mess of arguments, we take a max size and a Config struct to
	*   simulate named ctor parameters.  The Config struct has sensible defaults
	*   for everything, but is overloaded - if you specify a positive capacity,
	*   that will be used directly instead of computing it based on
	*   maxLoadFactor.
	*
	*   Create returns an LFHA::SmartPtr which is a unique_ptr with a custom
	*   deleter to make sure everything is cleaned up properly.
	*/
	struct Config
	{
		KeyT emptyKey{(KeyT)-1};
		KeyT lockedKey{(KeyT)-2};
		KeyT erasedKey{(KeyT)-3};
		double maxLoadFactor{0.8};
		double growthFactor{-1.0};
		uint32_t entryCountSize{1000};
		size_t capacity{0}; // if positive, overrides maxLoadFactor
		// Cannot have constexpr ctor because some compilers rightly complain.
		Config() = default;
	};

	//  Cannot have pre-instantiated const Config instance because of SIOF.
	static SmartPtr Create(size_t maxSize, const Config & c = Config())
	{
		GENERIC_ASSERT(1.0 > c.maxLoadFactor);
		GENERIC_ASSERT(0.0 < c.maxLoadFactor);
		GENERIC_ASSERT(c.emptyKey != c.lockedKey);
		size_t capacity = size_t(maxSize / c.maxLoadFactor);
		size_t sz = sizeof(LockFreeHashArray) + sizeof(value_type) * capacity;
		auto const mem = Allocator().allocate(sz);
		try {
			new (mem) LockFreeHashArray(
				capacity,
				c.emptyKey,
				c.lockedKey,
				c.erasedKey,
				c.maxLoadFactor,
				c.entryCountSize);
		} catch (...) {
			Allocator().deallocate(mem, sz);
			throw;
		}

		SmartPtr map(static_cast<LockFreeHashArray*>((void*)mem));

		/*
		* Mark all cells as empty.
		*
		* Note: we're bending the rules a little here accessing the key
		* element in our cells even though the cell object has not been
		* constructed, and casting them to atomic objects (see CellKeyPtr).
		* (Also, in fact we never actually invoke the value_type
		* constructor.)  This is in order to avoid needing to default
		* construct a bunch of value_type when we first start up: if you
		* have an expensive default constructor for the value type this can
		* noticeably speed construction time for an LFHA.
		*/
		for (size_t i = 0; i < map->capacity; ++i)
			CellKeyPtr(map->m_cells[i])->store(map->kEmptyKey, std::memory_order_relaxed);
		return map;
	}

	/*
	* find --
	*   Returns the iterator to the element if found, otherwise end().
	*   As an optional feature, the type of the key to look up (LookupKeyT) is
	*   allowed to be different from the type of keys actually stored (KeyT).
	*   This enables use cases where materializing the key is costly and usually
	*   redundant, e.g., canonicalizing/interning a set of strings and being able
	*   to look up by StringPiece. To use this feature, LookupHashFcn must take
	*   a LookupKeyT, and LookupEqualFcn must take KeyT and LookupKeyT as first
	*   and second parameter, respectively.
	*/
	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal>
	iterator find(const LookupKeyT & key)
	{
		return iterator(this, FindInternal<LookupKeyT, LookupHashFcn, LookupEqualFcn>(key).idx);
	}

	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal>
	const_iterator find(const LookupKeyT & k) const
	{
		return const_cast<LockFreeHashArray *>(this)->find<LookupKeyT, LookupHashFcn, LookupEqualFcn>(k);
	}

	/*
	* insert --
	*   Returns a pair with iterator to the element at r.first and bool success.
	*   Retrieve the index with ret.first.getIndex().
	*   Fails on key collision (does not overwrite) or if map becomes
	*   full, at which point no element is inserted, iterator is set to end(),
	*   and success is set false.  On collisions, success is set false, but the
	*   iterator is set to the existing entry.
	*/
	std::pair<iterator, bool> insert(const value_type & r)
	{
		return emplace(r.first, r.second);
	}

	std::pair<iterator, bool> insert(value_type && r)
	{
		return emplace(r.first, std::move(r.second));
	}

	/*
	* emplace --
	*   Same contract as insert(), but performs in-place construction
	*   of the value type using the specified arguments.
	*   Also, like find(), this method optionally allows 'key_in' to have a type
	*   different from that stored in the table; see find(). If and only if no
	*   equal key is already present, this method converts 'key_in' to a key of
	*   type KeyT using the provided LookupKeyToKeyFcn.
	*/
	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal,
		typename LookupKeyToKeyFcn = key_convert,
		typename... Args>
	std::pair<iterator, bool> emplace(LookupKeyT && keyIn, Args &&... args)
	{
		SimpleRetT ret = InsertInternal<
			LookupKeyT,
			LookupHashFcn,
			LookupEqualFcn,
			LookupKeyToKeyFcn>(keyIn, std::forward<Args>(args)...);
		return std::make_pair(iterator(this, ret.idx), ret.success);
	}

	// returns the number of elements erased - should never exceed 1
	size_t erase(KeyT keyIn)
	{
		GENERIC_ASSERT(keyIn != kEmptyKey);
		GENERIC_ASSERT(keyIn != kLockedKey);
		GENERIC_ASSERT(keyIn != kErasedKey);

		for (size_t idx = KeyToAnchorIdx(keyIn), numProbes = 0;;
			idx = ProbeFcn()(idx, numProbes, capacity)) {
			GENERIC_ASSERT(idx < capacity);
			value_type * cell = &m_cells[idx];
			KeyT currentKey = AcquireLoadKey(*cell);
			if (currentKey == kEmptyKey || currentKey == kLockedKey) {
				// If we hit an empty (or locked) element, this key does not exist. This
				// is similar to how it's handled in find().
				return 0;
			}
			if (EqualFcn()(currentKey, keyIn)) {
				// Found an existing entry for our key, attempt to mark it erased.
				// Some other thread may have erased our key, but this is ok.
				KeyT expect = currentKey;
				if (CellKeyPtr(*cell)->compare_exchange_strong(expect, kErasedKey)) {
					m_numErases.fetch_add(1, std::memory_order_relaxed);
					// Even if there's a value in the cell, we won't delete (or even
					// default construct) it because some other thread may be accessing it.
					// Locking it meanwhile won't work either since another thread may be
					// holding a pointer to it.
					// We found the key and successfully erased it.
					return 1;
				}
				// If another thread succeeds in erasing our key, we'll stop our search.
				return 0;
			}

			// NOTE: the way we count numProbes must be same in find(), insert(),
			// and erase(). Otherwise it may break probing.
			++numProbes;
			if (GENERIC_UNLIKELY(numProbes >= capacity))
				return 0;// probed every cell...faildestroy
		}
	}
	
	// clears all keys and values in the map and resets all counters.  Not thread safe.
	void clear()
	{
		for (size_t i = 0; i < capacity; ++i) {
			if (m_cells[i].first != kEmptyKey) {
				m_cells[i].~value_type();
				*const_cast<KeyT*>(&m_cells[i].first) = kEmptyKey;
			}
			GENERIC_ASSERT(m_cells[i].first == kEmptyKey);
		}
		m_numEntries.store(0);
		m_numPendingEntries.store(0);
		m_numErases.store(0, std::memory_order_relaxed);
		m_isFull.store(0, std::memory_order_relaxed);
	}

	size_t size() const
	{
		return m_numEntries.load() - m_numErases.load(std::memory_order_relaxed);
	}
	
  	bool empty() const { return size() == 0; }

	iterator begin()
	{
		iterator it(this, 0);
		it.AdvancePastEmpty();
		return it;
	}

	const_iterator begin() const
	{
		const_iterator it(this, 0);
		it.AdvancePastEmpty();
		return it;
	}

	iterator end() { return iterator(this, capacity); }

	const_iterator end() const { return const_iterator(this, capacity); }

	// See LockFreeHashMap::FindAt - access elements directly
	// WARNING: The following 2 functions will fail silently for hashtable
	// with capacity > 2^32
	iterator FindAt(uint32_t idx)
	{
		GENERIC_ASSERT(idx < capacity);
		return iterator(this, idx);
	}

	const_iterator FindAt(uint32_t idx) const
	{
		return const_cast<LockFreeHashArray*>(this)->FindAt(idx);
	}

	iterator makeIter(size_t idx) { return iterator(this, idx); }

	const_iterator makeIter(size_t idx) const { return const_iterator(this, idx); }

  	// The max load factor allowed for this map
	double MaxLoadFactor() const { return ((double)maxEntries) / capacity; }

private:
	inline void UnlockCell(const value_type * cell, KeyT newKey)
	{
		CellKeyPtr(*cell)->store(newKey, std::memory_order_release);
	}

	inline bool TryLockCell(const value_type * cell)
	{
		KeyT expect = kEmptyKey;
		return CellKeyPtr(*cell)->compare_exchange_strong(expect, kLockedKey, std::memory_order_acq_rel);
	}
  
	template <class LookupKeyT = key_type, class LookupHashFcn = hasher>
	inline size_t KeyToAnchorIdx(const LookupKeyT k) const
	{
		const size_t hashVal = LookupHashFcn()(k);
		const size_t probe = hashVal & m_kAnchorMask;
		return GENERIC_LIKELY(probe < capacity) ? probe : hashVal % capacity;
	}
  
	static KeyT RelaxedLoadKey(const value_type & r)
	{
		return CellKeyPtr(r)->load(std::memory_order_relaxed);
	}

	static KeyT AcquireLoadKey(const value_type & r)
	{
		return CellKeyPtr(r)->load(std::memory_order_acquire);
	}

	template <typename MaybeKeyT>
	void CheckLegalKeyIfKey(const MaybeKeyT & key)
	{
		detail::CheckLegalKeyIfKeyTImpl(key, kEmptyKey, kLockedKey, kErasedKey);
	}

	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal>
	SimpleRetT FindInternal(const LookupKeyT & keyIn)
	{
		ProbeFcn probeFcn;
		LookupEqualFcn leFcn;
		CheckLegalKeyIfKey<LookupKeyT>(keyIn);
		for (size_t idx = KeyToAnchorIdx<LookupKeyT, LookupHashFcn>(keyIn),
			 numProbes = 0;; idx = probeFcn(idx, numProbes, capacity)) {
			const KeyT key = AcquireLoadKey(m_cells[idx]);
			if (GENERIC_LIKELY(leFcn(key, keyIn)))
				return SimpleRetT(idx, true);
			if (GENERIC_UNLIKELY(key == kEmptyKey))
				return SimpleRetT(capacity, false);//hit an empty element, this key does not exist

			// NOTE: the way we count numProbes must be same in find(), insert(),
			// and erase(). Otherwise it may break probing.
			++numProbes;
			if (GENERIC_UNLIKELY(numProbes >= capacity))
				return SimpleRetT(capacity, false);// probed every cell...fail
		}
	}

  	template <
      typename LookupKeyT = key_type,
      typename LookupHashFcn = hasher,
      typename LookupEqualFcn = key_equal,
      typename LookupKeyToKeyFcn = detail::Identity,
      typename... Args>
	SimpleRetT InsertInternal(LookupKeyT keyIn, Args &&... args)
	{
		constexpr short NO_NEW_INSERTS = 1;
		constexpr short NO_PENDING_INSERTS = 2;
		CheckLegalKeyIfKey<LookupKeyT>(keyIn);

		size_t idx = KeyToAnchorIdx<LookupKeyT, LookupHashFcn>(keyIn);
		size_t numProbes = 0;
		for (;;) {
			GENERIC_ASSERT(idx < capacity);
			value_type * cell = &m_cells[idx];
			if (RelaxedLoadKey(*cell) == kEmptyKey) {
				// NOTE: m_isFull is set based on m_numEntries, so it's
				// possible to insert more than maxEntries entries. However, it's not
				// possible to insert past capacity.
				++m_numPendingEntries;
				if (m_isFull.load(std::memory_order_acquire)) {
					--m_numPendingEntries;

					// Before deciding whether this insert succeeded, this thread needs to
					// wait until no other thread can add a new entry.

					// Correctness assumes m_isFull is true at this point. If
					// another thread now does ++m_numPendingEntries, we expect it
					// to pass the m_isFull.load() test above. (It shouldn't insert
					// a new entry.)
					SpinWait([&] { return (m_isFull.load(std::memory_order_acquire) != NO_PENDING_INSERTS) && (m_numPendingEntries.load() != 0); });
					m_isFull.store(NO_PENDING_INSERTS, std::memory_order_release);

					// Don't insert past max load factor
					if (RelaxedLoadKey(*cell) == kEmptyKey)
						return SimpleRetT(capacity, false);
				} 
				else {
					// An unallocated cell. Try once to lock it. If we succeed, insert here.
					// If we fail, fall through to comparison below; maybe the insert that
					// just beat us was for this very key....
					if (TryLockCell(cell)) {
						KeyT keyNew;
						// Write the value - done before unlocking
						try {
							keyNew = LookupKeyToKeyFcn()(keyIn);
							using LookupKeyTNoConst = std::remove_const_t<LookupKeyT>;
							if constexpr (not std::is_same_v<KeyT, LookupKeyTNoConst>)
								CheckLegalKeyIfKey(keyNew);
							GENERIC_ASSERT(RelaxedLoadKey(*cell) == kLockedKey);
							// A const mapped_type is only constant once constructed, so cast
							// away any const for the placement new here.
							using mapped = typename std::remove_const_t<mapped_type>;
							new (const_cast<mapped*>(&cell->second)) ValueT(std::forward<Args>(args)...);
							UnlockCell(cell, keyNew); // Sets the new key
						} catch (...) {
							// Transition back to empty key---requires handling
							// locked->empty below.
							UnlockCell(cell, kEmptyKey);
							--m_numPendingEntries;
							throw;
						}
						// An erase() can race here and delete right after our insertion
						// Direct comparison rather than EqualFcn ok here. (we just inserted it)
						GENERIC_ASSERT(RelaxedLoadKey(*cell) == keyNew || RelaxedLoadKey(*cell) == kErasedKey);
						--m_numPendingEntries;
						++m_numEntries;
						if (m_numEntries.load() >= maxEntries)
							m_isFull.store(NO_NEW_INSERTS, std::memory_order_relaxed);
						return SimpleRetT(idx, true);
					}
					--m_numPendingEntries;
				}
			}
			GENERIC_ASSERT(RelaxedLoadKey(*cell) != kEmptyKey);
			if (kLockedKey == AcquireLoadKey(*cell))
				SpinWait([&]{ return kLockedKey == AcquireLoadKey(*cell); });

			const KeyT thisKey = AcquireLoadKey(*cell);
			if (LookupEqualFcn()(thisKey, keyIn)) {
				// Found an existing entry for our key, but we don't overwrite the previous value.
				return SimpleRetT(idx, false);
			} 
			else if (thisKey == kEmptyKey || thisKey == kLockedKey) {
				// We need to try again (i.e., don't increment numProbes or
				// advance idx): this case can happen if the constructor for
				// ValueT threw for this very cell (the rethrow block above).
				continue;
			}

			// NOTE: the way we count numProbes must be same in find(),
			// insert(), and erase(). Otherwise it may break probing.
			++numProbes;
			if (GENERIC_UNLIKELY(numProbes >= capacity))
				return SimpleRetT(capacity, false);// probed every cell...fail

			idx = ProbeFcn()(idx, numProbes, capacity);
		}
	}

private:
	std::atomic<size_t> m_numEntries{0};// Successful key inserts
	std::atomic<size_t> m_numPendingEntries{0}; // Used by insertInternal
	std::atomic<size_t> m_numErases{0}; // Successful key erases
	std::atomic<int64_t> m_isFull{0}; // Used by insertInternal
	value_type m_cells[0]; // This must be the last field of this class
};

/*
* LockFreeHashMap provides an interface somewhat similar to the
* UnorderedAssociativeContainer concept in C++.  This does not
* exactly match this concept (or even the basic Container concept),
* because of some restrictions imposed by our datastructure.
*
* Specific differences (there are quite a few):

* - Efficiently thread safe for inserts (main point of this stuff),
*   wait-free for lookups.
*
* - You can erase from this container, but the cell containing the key will
*   not be free or reclaimed.
*
* - You can erase everything by calling clear() (and you must guarantee only
*   one thread can be using the container to do that).
*
* - We aren't DefaultConstructible, CopyConstructible, Assignable, or
*   EqualityComparable.  (Most of these are probably not something
*   you actually want to do with this anyway.)
*
* - We don't support the various bucket functions, rehash(),
*   reserve(), or equal_range().  Also no constructors taking
*   iterators, although this could change.
*
* - Several insertion functions, notably operator[], are not
*   implemented.  It is a little too easy to misuse these functions
*   with this container, where part of the point is that when an
*   insertion happens for a new key, it will atomically have the
*   desired value.
*
* - The map has no templated insert() taking an iterator range, but
*   we do provide an insert(key, value).  The latter seems more
*   frequently useful for this container (to avoid sprinkling
*   make_pair everywhere), and providing both can lead to some gross
*   template error messages.
*
* - The Allocator must not be stateful (a new instance will be spun up for
*   each allocation), and its allocate() method must take a raw number of
*   bytes.
*
* - KeyT must be a 32 bit or 64 bit atomic integer type, and you must
*   define special 'locked' and 'empty' key values in the ctor
*
* - We don't take the Hash function object as an instance in the
*   constructor.
*
*/

template <
    class KeyT,
    class ValueT,
    class HashFcn,
    class EqualFcn,
    class Allocator,
    class ProbeFcn,
    class KeyConvertFcn>
class LockFreeHashMap
{
	using SubMap = 	LockFreeHashArray<KeyT, ValueT, HashFcn, EqualFcn, Allocator, ProbeFcn, KeyConvertFcn>;
public:
	using key_type = KeyT;
	using mapped_type = ValueT;
	using value_type = std::pair<const KeyT, ValueT>;
	using hasher = HashFcn;
	using key_equal = EqualFcn;
	using key_convert = KeyConvertFcn;
	using pointer = value_type*;
	using reference = value_type &;
	using const_reference = const value_type &;
	using difference_type = std::ptrdiff_t ;
	using size_type = std::size_t ;
	using Config = typename SubMap::Config ;

	template <class ContT, class IterVal, class SubIt>
	struct Iterator
	{
		Iterator() = default;

	private:
		friend class LockFreeHashMap;
		friend class detail::IteratorFacade<Iterator, IterVal, std::forward_iterator_tag>;
		explicit Iterator(ContT* lfhm, uint32_t subMap, const SubIt & subIt)
		: m_lfhm(lfhm), m_subMap(subMap), m_subIt(subIt) {}

		void Increment()
		{
			GENERIC_ASSERT(!isEnd());
			++m_subIt;
			CheckAdvanceToNextSubmap();
		}

		bool Equal(const Iterator & other) const
		{
			if (m_lfhm != other.m_lfhm) return false;

			if (isEnd() || other.isEnd())
				return isEnd() == other.isEnd();

			return m_subMap == other.m_subMap && m_subIt == other.m_subIt;
		}

	 	IterVal & Dereference() const { return *m_subIt; }

		bool isEnd() const { return m_lfhm == nullptr; }

		void CheckAdvanceToNextSubmap()
		{
			if (isEnd()) return;

			SubMap * thisMap = m_lfhm->m_subMaps[m_subMap].load(std::memory_order_relaxed);
			while (m_subIt == thisMap->end()) {
				// This sub iterator is done, advance to next one
				if (m_subMap + 1 < m_lfhm->m_numMapsAllocated.load(std::memory_order_acquire)) {
					++m_subMap;
					thisMap = m_lfhm->m_subMaps[m_subMap].load(std::memory_order_relaxed);
					m_subIt = thisMap->begin();
				} else {
					m_lfhm = nullptr;
					return;
				}
			}
		}
	private:
		ContT* m_lfhm{nullptr};
		uint32_t m_subMap{0};
		SubIt m_subIt{};
	};

	using iterator = Iterator<LockFreeHashMap, value_type, typename SubMap::iterator>;
	using const_iterator = Iterator<const LockFreeHashMap, const value_type, typename SubMap::const_iterator>;

public:
 	const float kGrowthFrac{0}; // How much to grow when we run out of capacity.

	// The constructor takes a finalSizeEst which is the optimal
	// number of elements to maximize space utilization and performance,
	// and a Config object to specify more advanced options.
  	explicit LockFreeHashMap(size_t finalSizeEst, const Config & c = Config())
     : kGrowthFrac(c.growthFactor < 0 ? 1.0f - c.maxLoadFactor : c.growthFactor)
	{
		GENERIC_ASSERT(c.maxLoadFactor > .0f && c.maxLoadFactor < 1.0f);
		m_subMaps[0].store(SubMap::Create(finalSizeEst, c).release(), std::memory_order_relaxed);
		auto subMapCount = kNumSubMaps;
		for (size_t i = 1; i < subMapCount; ++i)
			m_subMaps[i].store(nullptr, std::memory_order_relaxed);
		m_numMapsAllocated.store(1, std::memory_order_relaxed);
	}

	LockFreeHashMap(const LockFreeHashMap &) = delete;
	LockFreeHashMap& operator=(const LockFreeHashMap &) = delete;

	~LockFreeHashMap()
	{
		const unsigned int numMaps = m_numMapsAllocated.load(std::memory_order_relaxed);
		for(size_t i = 0; i < numMaps; ++i) {
			SubMap * thisMap = m_subMaps[i].load(std::memory_order_relaxed);
			GENERIC_ASSERT(thisMap);
			SubMap::Destroy(thisMap);
		}
	}

	key_equal key_eq() const { return key_equal(); }

	hasher hash_function() const { return hasher(); }

	/*
	* insert --
	*
	*   Returns a pair with iterator to the element at r.first and
	*   success.  Retrieve the index with ret.first.getIndex().
	*
	*   Does not overwrite on key collision, but returns an iterator to
	*   the existing element (since this could due to a race with
	*   another thread, it is often important to check this return
	*   value).
	*
	*   Allocates new sub maps as the existing ones become full.  If
	*   all sub maps are full, no element is inserted, and
	*   full error is thrown.
	*/
	std::pair<iterator, bool> insert(const value_type & r)
	{
		return emplace(r.first, r.second);
	}
	
	std::pair<iterator, bool> insert(key_type k, const mapped_type & v)
	{
		return emplace(k, v);
	}

	std::pair<iterator, bool> insert(value_type && r)
	{
		return emplace(r.first, std::move(r.second));
	}

	std::pair<iterator, bool> insert(key_type k, mapped_type && v)
	{
		return emplace(k, std::move(v));
	}

	/*
	* emplace --
	*
	*   Same contract as insert(), but performs in-place construction
	*   of the value type using the specified arguments.
	*
	*   Also, like find(), this method optionally allows 'keyIn' to have a type
	*   different from that stored in the table; see find(). If and only if no
	*   equal key is already present, this method converts 'keyIn' to a key of
	*   type KeyT using the provided LookupKeyToKeyFcn.
	*/
	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal,
		typename LookupKeyToKeyFcn = key_convert,
		typename... Args>
	std::pair<iterator, bool> emplace(LookupKeyT k, Args &&... args)
	{
		auto ret = InsertInternal<
			LookupKeyT,
			LookupHashFcn,
			LookupEqualFcn,
			LookupKeyToKeyFcn>(k, std::forward<Args>(args)...);
		SubMap* subMap = m_subMaps[ret.i].load(std::memory_order_relaxed);
		return std::make_pair(iterator(this, ret.i, subMap->makeIter(ret.j)), ret.success);
	}

	/*
	* find --
	*
	*   Returns the iterator to the element if found, otherwise end().
	*
	*   As an optional feature, the type of the key to look up (LookupKeyT) is
	*   allowed to be different from the type of keys actually stored (KeyT).
	*
	*   This enables use cases where materializing the key is costly and usually
	*   redundant, e.g., canonicalizing/interning a set of strings and being able
	*   to look up by StringPiece. To use this feature, LookupHashFcn must take
	*   a LookupKeyT, and LookupEqualFcn must take KeyT and LookupKeyT as first
	*   and second parameter, respectively.
	*/
	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal>
	iterator find(const LookupKeyT & k)
	{
		SimpleRetT ret = FindInternal<LookupKeyT, LookupHashFcn, LookupEqualFcn>(k);
		if (not ret.success) return end();
		
		SubMap * subMap = m_subMaps[ret.i].load(std::memory_order_relaxed);
		return iterator(this, ret.i, subMap->makeIter(ret.j));
	}

	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal>
	const_iterator find(const LookupKeyT & k) const
	{
		return const_cast<LockFreeHashMap*>(this)->find<LookupKeyT, LookupHashFcn, LookupEqualFcn>(k);
	}

	/*
	* findAt --
	*
	*   Returns an iterator into the map.
	*
	*   idx should only be an unmodified value returned by calling getIndex() on
	*   a valid iterator returned by find() or insert(). If idx is invalid you
	*   have a bug and the process aborts.
	*/
	iterator FindAt(uint32_t idx)
	{
		SimpleRetT ret = FindAtInternal(idx);
		GENERIC_ASSERT(ret.i < numSubMaps());
		return iterator(this, ret.i, m_subMaps[ret.i].load(std::memory_order_relaxed)->makeIter(ret.j));
	}

	const_iterator FindAt(uint32_t idx) const
	{
		return const_cast<LockFreeHashMap*>(this)->FindAt(idx);
	}

	// Number of sub maps allocated so far to implement this map.  The more there
	// are, the worse the performance.
	int numSubMaps() const
	{
		return m_numMapsAllocated.load(std::memory_order_acquire);
	}

	iterator begin()
	{
		iterator it(this, 0, m_subMaps[0].load(std::memory_order_relaxed)->begin());
		it.CheckAdvanceToNextSubmap();
		return it;
	}

	const_iterator begin() const
	{
		const_iterator it(this, 0, m_subMaps[0].load(std::memory_order_relaxed)->begin());
		it.CheckAdvanceToNextSubmap();
		return it;
	}

	iterator end() { return iterator(); }

	const_iterator end() const { return const_iterator(); }

	/*
	* erase --
	*   Erases key k from the map
	*   Returns 1 if the key is found and erased, and 0 otherwise.
	*/
	size_type erase(key_type k)
	{
		int const numMaps = m_numMapsAllocated.load(std::memory_order_acquire);
		for (size_t i = 0; i < numMaps; ++i) {
			// Check each map successively.  If one succeeds, we're done!
			if (m_subMaps[i].load(std::memory_order_relaxed)->erase(k)) return 1;
		}
		// Didn't find our key...
		return 0;
	}

	/*
	* size --
	*
	*  Returns the exact size of the map.  Note this is not as cheap as typical
	*  size() implementations because, for each LockFreeHashArray in this LFAM,
	*  we need to grab a lock and accumulate the values from all the thread local
	*  counters.
	*/
	size_t size() const
	{
		size_t totalSize(0);
		int const numMaps = m_numMapsAllocated.load(std::memory_order_acquire);
		for (size_t i = 0; i < numMaps; ++i)
			totalSize += m_subMaps[i].load(std::memory_order_relaxed)->size();
		return totalSize;
	}

  	bool empty() const { return size() == 0; }

  	size_type count(const key_type & k) const { return find(k) == end() ? 0 : 1; }

	// Total capacity - summation of capacities of all submaps.
	size_t capacity() const
	{
		size_t totalCap(0);
		int const numMaps = m_numMapsAllocated.load(std::memory_order_acquire);
		for (size_t i = 0; i < numMaps; ++i) {
			totalCap += m_subMaps[i].load(std::memory_order_relaxed)->capacity;
		}
		return totalCap;
	}

	// Number of new insertions until current submaps are all at max load factor.
  	size_t SpaceRemaining() const
	{
		size_t spaceRem(0);
		int const numMaps = m_numMapsAllocated.load(std::memory_order_acquire);
		for (size_t i = 0; i < numMaps; ++i) {
			SubMap * thisMap = m_subMaps[i].load(std::memory_order_relaxed);
			spaceRem += std::max(0, thisMap->maxEntries - &thisMap->m_numEntries.load(std::memory_order_acquire));
		}
		return spaceRem;
	}

	/*
	* clear --
	*
	*   Wipes all keys and values from primary map and destroys all secondary
	*   maps.  Primary map remains allocated and thus the memory can be reused
	*   in place.  Not thread safe.
	*/
	void clear()
	{
		m_subMaps[0].load(std::memory_order_relaxed)->clear();
		int const numMaps = m_numMapsAllocated.load(std::memory_order_relaxed);
		for (size_t i = 1; i < numMaps; ++i) {
			SubMap* thisMap = m_subMaps[i].load(std::memory_order_relaxed);
			GENERIC_ASSERT(thisMap);
			SubMap::Destroy(thisMap);
			m_subMaps[i].store(nullptr, std::memory_order_relaxed);
		}
		m_numMapsAllocated.store(1, std::memory_order_relaxed);
	}

	/* Advanced functions for direct access: */
	uint32_t Rec2Idx(const value_type & r, bool mayInsert = true)
	{
		SimpleRetT ret = mayInsert ? InsertInternal(r.first, r.second) : FindInternal(r.first);
		return EncodeIndex(ret.i, ret.j);
	}

	uint32_t Rec2Idx(value_type && r, bool mayInsert = true)
	{
		SimpleRetT ret = mayInsert ? InsertInternal(r.first, std::move(r.second)) : FindInternal(r.first);
		return EncodeIndex(ret.i, ret.j);
	}

	uint32_t Rec2Idx(key_type k, const mapped_type & v, bool mayInsert = true)
	{
		SimpleRetT ret = mayInsert ? InsertInternal(k, v) : FindInternal(k);
		return EncodeIndex(ret.i, ret.j);
	}

	uint32_t Rec2Idx(key_type k, mapped_type && v, bool mayInsert = true) 
	{
		SimpleRetT ret = mayInsert ? InsertInternal(k, std::move(v)) : FindInternal(k);
		return EncodeIndex(ret.i, ret.j);
	}

	uint32_t Key2Idx(const KeyT k, bool mayInsert = false)
	{
		return Rec2Idx(value_type(k), mayInsert);
	}

	const value_type & Idx2Rec(uint32_t idx) const
	{
		SimpleRetT ret = FindAtInternal(idx);
		return m_subMaps[ret.i].load(std::memory_order_relaxed)->Idx2Rec(ret.j);
	}

private:
	// This limits primary submap size to 2^31 ~= 2 billion, secondary submap
	// size to 2^(32 - kNumSubMapBits_ - 1) = 2^27 ~= 130 million, and num subMaps
	// to 2^kNumSubMapBits = 16.
	static constexpr uint32_t kNumSubMapBits = 4;
	static constexpr uint32_t kSecondaryMapBit = 1u << 31; // Highest bit
	static constexpr uint32_t kSubMapIndexShift = 32 - kNumSubMapBits - 1;
	static constexpr uint32_t kSubMapIndexMask = (1 << kSubMapIndexShift) - 1;
	static constexpr uint32_t kNumSubMaps = 1 << kNumSubMapBits;
	static constexpr uintptr_t kLockedPtr = 0x88ULL << 48; // invalid pointer

	struct SimpleRetT
	{
		uint32_t i;
		size_t j;
		bool success;
		SimpleRetT(uint32_t ii, size_t jj, bool s) : i(ii), j(jj), success(s) {}
		SimpleRetT() = default;
	};

	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal,
		typename LookupKeyToKeyFcn = key_convert,
		typename... Args>
	SimpleRetT InsertInternal(LookupKeyT key, Args &&... args)
	{
		while (true) {
			auto nextMapIdx = m_numMapsAllocated.load(std::memory_order_acquire);// this maintains our state
			typename SubMap::SimpleRetT ret;
			for (size_t i = 0; i < nextMapIdx; ++i) {
				// insert in each map successively.  If one succeeds, we're done!
				SubMap* subMap = m_subMaps[i].load(std::memory_order_relaxed);
				ret = subMap->template InsertInternal<
					LookupKeyT,
					LookupHashFcn,
					LookupEqualFcn,
					LookupKeyToKeyFcn>(key, std::forward<Args>(args)...);
				if (ret.idx == subMap->capacity_)
					continue; // map is full, so try the next one
				// Either collision or success - insert in either case
				return SimpleRetT(i, ret.idx, ret.success);
			}

			// If we made it this far, all maps are full and we need to try to allocate
			// the next one.
			SubMap* primarySubMap = m_subMaps[0].load(std::memory_order_relaxed);
			if (nextMapIdx >= kNumSubMaps || primarySubMap->capacity * kGrowthFrac < 1.0) {
				// Can't allocate any more sub maps.
				ThrowException("hash map is full");
			}

			if (TryLockMap(nextMapIdx)) {
				// Alloc a new map and shove it in.  We can change whatever
				// we want because other threads are waiting on us...
				size_t numCellsAllocated = (size_t)(primarySubMap->capacity * std::pow(1.0 + kGrowthFrac, nextMapIdx - 1));
				size_t newSize = size_t(numCellsAllocated * kGrowthFrac);
				GENERIC_ASSERT(m_subMaps[nextMapIdx].load(std::memory_order_relaxed) == (SubMap*)kLockedPtr);
				// create a new map using the settings stored in the first map

				Config config;
				config.emptyKey = primarySubMap->kEmptyKey;
				config.lockedKey = primarySubMap->kLockedKey;
				config.erasedKey = primarySubMap->kErasedKey;
				config.maxLoadFactor = primarySubMap->maxLoadFactor();
				config.entryCountSize = primarySubMap->m_numEntries.load(std::memory_order_relaxed);
				m_subMaps[nextMapIdx].store(SubMap::Create(newSize, config).release(), std::memory_order_relaxed);

				// Publish the map to other threads.
				m_numMapsAllocated.fetch_add(1, std::memory_order_release);
				GENERIC_ASSERT(nextMapIdx + 1 == m_numMapsAllocated.load(std::memory_order_relaxed));
			} else {
				// If we lost the race, we'll have to wait for the next map to get
				// allocated before doing any insertion here.
				SpinWait([&] { return nextMapIdx >= m_numMapsAllocated.load(std::memory_order_acquire); });
			}

			// Relaxed is ok here because either we just created this map, or we
			// just did a spin wait with an acquire load on m_numMapsAllocated.
			SubMap * loadedMap = m_subMaps[nextMapIdx].load(std::memory_order_relaxed);
			GENERIC_ASSERT(loadedMap && loadedMap != (SubMap *)kLockedPtr);
			ret = loadedMap->InsertInternal(key, std::forward<Args>(args)...);
			if (ret.idx != loadedMap->capacity)
				return SimpleRetT(nextMapIdx, ret.idx, ret.success);
	
			// We took way too long and the new map is already full...try again from
			// the top (this should pretty much never happen).
		}
	}

	template <
		typename LookupKeyT = key_type,
		typename LookupHashFcn = hasher,
		typename LookupEqualFcn = key_equal>
	SimpleRetT FindInternal(const LookupKeyT k) const
	{
		const SubMap * primaryMap = m_subMaps[0].load(std::memory_order_relaxed);
		auto ret = primaryMap->template FindInternal<LookupKeyT, LookupHashFcn, LookupEqualFcn>(k);
		if (GENERIC_LIKELY(ret.idx != primaryMap->capacity))
			return SimpleRetT(0, ret.idx, ret.success);
		
		const unsigned int numMaps = m_numMapsAllocated.load(std::memory_order_acquire);
		for (size_t i = 1; i < numMaps; ++i) {
			// Check each map successively.  If one succeeds, we're done!
			SubMap * thisMap = m_subMaps[i].load(std::memory_order_relaxed);
			ret = thisMap->template findInternal<LookupKeyT, LookupHashFcn, LookupEqualFcn>(k);
			if (GENERIC_LIKELY(ret.idx != thisMap->capacity_))
				return SimpleRetT(i, ret.idx, ret.success);
		}
		// Didn't find our key...
		return SimpleRetT(numMaps, 0, false);
	}

	SimpleRetT FindAtInternal(uint32_t idx) const
	{
		uint32_t subMapIdx, subMapOffset;
		if (idx & kSecondaryMapBit) {
			// idx falls in a secondary map
			idx &= ~kSecondaryMapBit; // unset secondary bit
			subMapIdx = idx >> kSubMapIndexShift;
			GENERIC_ASSERT(subMapIdx < m_numMapsAllocated.load(std::memory_order_relaxed));
			subMapOffset = idx & kSubMapIndexMask;
		} else {
			// idx falls in primary map
			subMapIdx = 0;
			subMapOffset = idx;
		}
		return SimpleRetT(subMapIdx, subMapOffset, true);
	}

	bool TryLockMap(unsigned int idx)
	{
		SubMap * val = nullptr;
		return m_subMaps[idx].compare_exchange_strong(val, (SubMap*)kLockedPtr, std::memory_order_acquire);
	}

	// EncodeIndex -- Encode the submap index and offset into return.
	// index_ret must be pre-populated with the submap offset.
	//
	// We leave index_ret untouched when referring to the primary map
	// so it can be as large as possible (31 data bits).  Max size of
	// secondary maps is limited by what can fit in the low 27 bits.
	//
	// Returns the following bit-encoded data in index_ret:
	//   if subMap == 0 (primary map) =>
	//     bit(s)          value
	//         31              0
	//       0-30  submap offset (index_ret input)
	//
	//   if subMap > 0 (secondary maps) =>
	//     bit(s)          value
	//         31              1
	//      27-30   which subMap
	//       0-26  subMap offset (index_ret input)
  	static uint32_t EncodeIndex(uint32_t subMap, uint32_t offset)
	{
		GENERIC_ASSERT((offset & kSecondaryMapBit) == 0); // offset can't be too big
		if (subMap == 0) return offset;

		// Make sure subMap isn't too big
		GENERIC_ASSERT((subMap >> kNumSubMapBits) == 0);
		// Make sure subMap bits of offset are clear
		GENERIC_ASSERT((offset & (~kSubMapIndexMask | kSecondaryMapBit)) == 0);

		// Set high-order bits to encode which submap this index belongs to
		return offset | (subMap << kSubMapIndexShift) | kSecondaryMapBit;
	}

private:
  std::atomic<SubMap*> m_subMaps[kNumSubMaps];
  std::atomic<uint32_t> m_numMapsAllocated;
};

template <
    class KeyT,
    class ValueT,
    class HashFcn = std::hash<KeyT>,
    class EqualFcn = std::equal_to<KeyT>,
    class Allocator = std::allocator<char>>
using QuadraticProbLockFreeHashMap = LockFreeHashMap<KeyT, ValueT, HashFcn, EqualFcn, Allocator, LockFreeHashArrayQuadraticProbeFcn>;

} // namespace generic::thread