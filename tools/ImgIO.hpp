/**
 * @file ImgIO.hpp
 * @author bwu
 * @brief image i/o
 * @version 0.1
 * @date 2023-08-28
 */
#pragma once
#include <string_view>
#include <type_traits>
#include <vector>
namespace generic::img {

namespace png {

inline static bool save(std::string_view filename, size_t w, size_t h, unsigned char * img, int alpha);

namespace detail {

class DataChunk
{
public:
/**
 * @brief Brief description of DataChunk.
 * @param startSize
 * @return explicit
 */
    explicit DataChunk(size_t startSize) { m_bytes.reserve(startSize); }

/**
 * @brief Brief description of Add.
 * @param str
 * @return void
 */
    void Add(std::string_view str) { for (auto c : str) m_bytes.emplace_back(c); }

    template <typename... Chars>
/**
 * @brief Brief description of AddChars.
 * @param chars
 * @return void
 */
    void AddChars(Chars&&... chars) { (m_bytes.emplace_back(chars), ...); }

    template <typename T, typename = std::enable_if_t<std::is_integral_v<T> > >
    void Add(T val) = delete;

/**
 * @brief Brief description of Add.
 * @param val
 * @return void
 */
    void Add(std::byte val) { m_bytes.emplace_back(static_cast<unsigned char>(val)); }

/**
 * @brief Brief description of Add.
 * @param val
 * @return void
 */
    void Add(uint16_t val) { AddChars(val, val >> 8); }

/**
 * @brief Brief description of Add.
 * @param val
 * @return void
 */
    void Add(uint32_t val) { AddChars(val >> 24, val >> 16, val >> 8, val); }

    void AddCRC()
    {
        const auto kIgnoredBytesCount{4};
        Add(calculateCRC(m_bytes.begin() + kIgnoredBytesCount, m_bytes.end()));
    }

    void UpdateAddler(unsigned char u)
    {
        m_a = (m_a + (u)) % 65521;
        m_b = (m_b + m_a) % 65521;
    }

    void AddAddler()
    {
        uint32_t addler{(m_b << 16) | m_a};
        Add(addler);
    }

/**
 * @brief Brief description of Data.
 * @return std::vector<unsigned char> &
 */
    std::vector<unsigned char> & Data() { return m_bytes; }
    const std::vector<unsigned char> & Data() const { return m_bytes; }

private:
    template <typename ForwardIt>
    uint32_t calculateCRC(ForwardIt first, ForwardIt last)
    {
        constexpr std::array<uint32_t, 16> kCrc{0x0,        0x1db71064, 0x3b6e20c8, 0x26d930ac,
                                                0x76dc4190, 0x6b6b51f4, 0x4db26158, 0x5005713c,
                                                0xedb88320, 0xf00f9344, 0xd6d6a3e8, 0xcb61b38c,
                                                0x9b64c2b0, 0x86d3d2d4, 0xa00ae278, 0xbdbdf21c};
        uint32_t crc = ~0U;
        while (first != last) {
            crc ^= (*first++);
            crc = (crc >> 4) ^ kCrc[crc & 15];
            crc = (crc >> 4) ^ kCrc[crc & 15];
        }
        return ~crc;
    }    
private:
    std::vector<unsigned char> m_bytes;
    uint32_t m_a{1};
    uint32_t m_b{0};
};

inline static std::vector<unsigned char> createIHDRChunk(uint32_t w, uint32_t h, int alpha)
{
    std::byte const kDpeth{8};
    uint32_t const kChunkSize{13};
    DataChunk hdrChunk(26);
    hdrChunk.Add(kChunkSize);
    hdrChunk.Add("IHDR");
    hdrChunk.Add(w);
    hdrChunk.Add(h);
    hdrChunk.Add(kDpeth);
    if (alpha == 1) hdrChunk.Add(std::byte(6));
    else hdrChunk.Add(std::byte(2));
    hdrChunk.AddChars('\0', '\0', '\0');
    hdrChunk.AddCRC();
    return hdrChunk.Data();
}

inline static std::vector<unsigned char> CreateDataChunk(uint32_t w, uint32_t h, int alpha, unsigned char * data)
{
    uint32_t const kPixelWidth = alpha ? 4 : 3;
    uint16_t const kRowWidth = w * kPixelWidth + 1;
    uint32_t kChunkSize = 2 + h * (5 + kRowWidth) + 4;

    DataChunk dataChunk(kChunkSize * 2);
    dataChunk.Add(kChunkSize);    // 4
    dataChunk.Add("IDAT\x78\1");  // 6 / 10

    for (uint32_t i = 0; i < h; ++i) {
        if (i != h - 1)
        dataChunk.Add(std::byte(0));  // 1 / 11
        else
        dataChunk.Add(std::byte(1));
        dataChunk.Add(kRowWidth);                          // 2 / 13
        dataChunk.Add(static_cast<uint16_t>(~kRowWidth));  // 2 / 15

        unsigned char const kNoFilter = 0;
        dataChunk.Add(std::byte(kNoFilter));
        dataChunk.UpdateAddler(kNoFilter);
        for (uint16_t i = 0; i < kRowWidth - 1; ++i, data++) {
            dataChunk.Add(std::byte(*data));
            dataChunk.UpdateAddler(*data);
        }
    }

    dataChunk.AddAddler();
    dataChunk.AddCRC();
    return dataChunk.Data();
}

} // namespace detail


} // namspace png

} // namespace generic::img
