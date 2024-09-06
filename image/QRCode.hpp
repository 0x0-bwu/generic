/**
 * @file Exception.hpp
 * @author bwu
 * @brief QRcode, modified from https://github.com/soleilpqd/QRMatrix.CPP
 * @version 0.1
 * @date 2024-08-30
 */
#pragma once
#include "generic/common/Exception.hpp"
#include "generic/tools/FileSystem.hpp"
#include "generic/common/System.hpp"

#include <boost/gil/extension/io/png.hpp>
#include <boost/gil.hpp>
#include <array>
#include <map>
namespace generic::img::qr
{

using UnsignedByte = unsigned char;
using Unsigned2Bytes = unsigned short;
using Unsigned4Bytes = size_t;

inline static constexpr int QR_MAX_VERSION = 40;
inline static constexpr int QR_VERSION_OFFSET = 4;
inline static constexpr int QR_MIN_DIMENSION = 21;
inline static constexpr int MICROQR_MAX_VERSION = 4;
inline static constexpr int MICROQR_VERSION_OFFSET = 2;
inline static constexpr int MICROQR_MIN_DIMENSION = 11;
inline static constexpr size_t DEFAULT_ECI_ASSIGMENT_VALUE = 3;
enum class EncodingMode
{
    NUMERIC = 0b0001,// 0-9
    ALPHA_NUMERIC = 0b0010,// 0-9, A-Z, $, %, *, +, -, ., /, :, ' '
    BYTE = 0b0100,// ISO-8859-1 (Unicode Latin-1)
    KANJI = 0b1000,// Shift JIS double-byte characters
};

enum class ErrorCorrectionLevel
{
    LOW = 0b01,//7%
    MEDIUM = 0b00,//15%
    QUARTER = 0b11,//25%
    HIGH = 0b10,//30%
};

namespace detail {

inline static constexpr size_t NUM_TRIPLE_DIGITS_BITS_LEN = 10;
inline static constexpr size_t NUM_DOUBLE_DIGITS_BITS_LEN = 7;
inline static constexpr size_t NUM_SINGLE_DIGIT_BITS_LEN = 4;
inline static constexpr size_t ALPHA_NUM_MULTIPLICATION = 45;
inline static constexpr size_t ALPHA_NUM_PAIR_CHARS_BITS_LEN = 11;
inline static constexpr size_t ALPHA_NUM_SINGLE_CHAR_BITS_LEN = 6;

enum BoardCell {
        /// Cell is not filled yet
        NEUTRAL          = 0x00,
        /// Cell is white (lower 4 bits)
        UNSET            = 0x05,
        /// Cell is black (lower 4 bits)
        SET              = 0x0A,
        /// Cell is finder (higher 4 bits)
        FINDER           = 0x10,
        /// Cell is separator (higher 4 bits)
        SEPARATOR        = 0x20,
        /// Cell is alignment (higher 4 bits)
        ALIGNMENT        = 0x30,
        /// Cell is timing (higher 4 bits)
        TIMING           = 0x40,
        /// Cell is dark module (higher 4 bits)
        DARK             = 0x50,
        /// Cell is format (higher 4 bits)
        FORMAT           = 0x60,
        /// Cell is version (higher 4 bits) (version ≥ 7)
        VERSION          = 0x70,
        /// Cell is ErrorCorrection
        ERROR_CORRECTION = 0x80,
        /// Cell is Remainder
        REMAINDER        = 0x90,
        /// cell & lowMask = lower 4 bits value
        LOWMASK          = 0x0F,
        /// cell & highMask = higher 4 bits value
        HIGHMASK         = 0xF0,
        /// cell & funcMask > 0 => cell is function module
        FUNCMASK         = 0x70,
};

inline int ConvertErrorCorrectionLevelToIndex(ErrorCorrectionLevel level) {
    switch (level) {
    case ErrorCorrectionLevel::LOW:
        return 0;
    case ErrorCorrectionLevel::MEDIUM:
        return 1;
    case ErrorCorrectionLevel::QUARTER:
        return 2;
    case ErrorCorrectionLevel::HIGH:
        return 3;
    }
    return -1;
}

 /// Get QR dimension by its version (version in 1...40 ~ dimension 21...177).
inline UnsignedByte GetDimensionByVersion(UnsignedByte version)
{
    return (version - 1) * QR_VERSION_OFFSET + QR_MIN_DIMENSION;
}

/// Get MicroQR dimension by its version (version in 1...4 ~ dimension 11...17).
inline UnsignedByte GetMicroDimensionByVersion(UnsignedByte version)
{
    return (version - 1) * MICROQR_VERSION_OFFSET + MICROQR_MIN_DIMENSION;
}

/// Number of bits for character counts indicator.
inline size_t GetCharactersCountIndicatorLength(UnsignedByte version, EncodingMode mode)
{
    if (1 <= version && version <= 9) {
        switch (mode) {
        case EncodingMode::NUMERIC:
            return 10;
        case EncodingMode::ALPHA_NUMERIC:
            return 9;
        case EncodingMode::BYTE:
            return 8;
        case EncodingMode::KANJI:
            return 8;
        default:
            return 0;
        }
    } 
    else if (10 <= version && version <= 26) {
         switch (mode) {
        case EncodingMode::NUMERIC:
            return 12;
        case EncodingMode::ALPHA_NUMERIC:
            return 11;
        case EncodingMode::BYTE:
            return 16;
        case EncodingMode::KANJI:
            return 10;
        default:
            return 0;
        }
    } 
    else if (27 <= version && version <= 40) {
        switch (mode) {
        case EncodingMode::NUMERIC:
            return 14;
        case EncodingMode::ALPHA_NUMERIC:
            return 13;
        case EncodingMode::BYTE:
            return 16;
        case EncodingMode::KANJI:
            return 12;
        default:
            return 0;
        }
    }
    return 0;
}

/// Number of bits of character counts indicator for MicroQR.
inline size_t GetMicroCharactersCountIndicatorLength(UnsignedByte version, EncodingMode mode)
{
    switch (mode) {
    case EncodingMode::NUMERIC:
        return version + 2;
    case EncodingMode::ALPHA_NUMERIC:
        return (version < 2) ? 0 : version + 1;
    case EncodingMode::BYTE:
        return (version < 3) ? 0 : version + 1;
    case EncodingMode::KANJI:
        return (version < 3) ? 0 : version;
    default:
        return 0;
    }
    return 0;
}

/// Number of bits of mode indicator for MicroQR.
inline size_t GetMicroModeIndicatorLength(UnsignedByte version, EncodingMode mode)
{
    if (version < 1 || version > MICROQR_MAX_VERSION) {
        return 0;
    }
    switch (mode) {
    case EncodingMode::NUMERIC:
    case EncodingMode::ALPHA_NUMERIC:
        return (version < 2) ? 0 : version - 1;
    case EncodingMode::BYTE:
    case EncodingMode::KANJI:
        return (version < 3) ? 0 : version - 1;
    default:
        return 0;
    } 
}

/// Number of bits of terminator for MicroQR.
inline size_t GetMicroTerminatorLength(UnsignedByte version)
{
    if (version < 1 || version > MICROQR_MAX_VERSION) {
        return 0;
    }
    return version * 2 + 1;
}

/// Map Encoding mode value to MicroQR
inline UnsignedByte GetMicroQREncodingModeValue(EncodingMode mode)
{
    switch (mode) {
    case EncodingMode::NUMERIC:
        return 0;
    case EncodingMode::ALPHA_NUMERIC:
        return 1;
    case EncodingMode::BYTE:
        return 2;
    case EncodingMode::KANJI:
        return 3;
    default:
        return 0;
    }
}

/// Map ErrorCorrectionLevel value to MicroQR
inline UnsignedByte GetMicroQRErrorCorrectionLevelValue(ErrorCorrectionLevel level, UnsignedByte version)
{
    switch (version) {
    case 1:
        return 0;
    case 2:
        switch (level) {
        case ErrorCorrectionLevel::LOW:
            return 0b001;
        case ErrorCorrectionLevel::MEDIUM:
            return 0b010;
        default:
            break;
        }
    case 3:
        switch (level) {
        case ErrorCorrectionLevel::LOW:
            return 0b011;
        case ErrorCorrectionLevel::MEDIUM:
            return 0b100;
        default:
            break;
        }
    case 4:
        switch (level) {
        case ErrorCorrectionLevel::LOW:
            return 0b101;
        case ErrorCorrectionLevel::MEDIUM:
            return 0b110;
        case ErrorCorrectionLevel::QUARTER:
            return 0b111;
        default:
            break;
        }
    default:
        break;
    }
    return 0;    
}

/// Allocate `count` of bytes (and reset to 0)
inline UnsignedByte * Allocate(size_t count) //todo remove after 
{
    UnsignedByte * buffer = new UnsignedByte[count];
    for (size_t idx = 0; idx < count; idx++) {
        buffer[idx] = 0;
    }
    return buffer;  
}

/// Copy `count` bits of source BYTE (from bit `sourceStartIndex`) into `destination` BYTE stating from bit at `destStartIndex`.
/// Bits index is from left to right.
/// Bit index must be 0...7 (because this copies 1 byte only).
/// `count` must be 1...8 (bits - size of `source`), also `count` must be < 8 - destStartIndex && 8 - sourceStartIndex.
inline void CopyBits(UnsignedByte sourceByte, size_t sourceStartIndex, UnsignedByte* destination, size_t destStartIndex, size_t count)
{
    if (sourceStartIndex >= 8 || destStartIndex >= 8) {
        ThrowException("bit index must be 0...7.");
    }
    if (count == 0 || count > (8 - destStartIndex) || count > (8 - sourceStartIndex)) {
        ThrowException("Count must be > 0 and < 8 - sourceStartIndex and < 8 - destStartIndex.");
    }
    // Pattern for bits from source
    UnsignedByte pattern = 0x80;
    for (size_t idx = 1; idx < count; idx++) {
        pattern = (pattern >> 1) | 0x80;
    }
    pattern = pattern >> sourceStartIndex;
    // Clean source
    UnsignedByte cleanSource = sourceByte & pattern;
    // Reposition
    if (sourceStartIndex > destStartIndex) {
        size_t offset = sourceStartIndex - destStartIndex;
        cleanSource = cleanSource << offset;
        pattern = pattern << offset;
    } 
    else if (sourceStartIndex < destStartIndex) {
        size_t offset = destStartIndex - sourceStartIndex;
        cleanSource = cleanSource >> offset;
        pattern = pattern >> offset;
    }
    // Clean destination
    pattern = ~pattern;
    *destination = *destination & pattern;
    // apply
    *destination = *destination | cleanSource;   
}

/// Copy `count` bits of source (from bit 0th) into `destination` stating from bit at `startIndex`.
inline void CopyBits(UnsignedByte* source, size_t sourceLength, size_t sourceStartIndex, bool isSourceOrderReversed, 
                            UnsignedByte* destination, size_t destStartIndex, size_t count)
{
    if (count == 0 || count > sourceLength * 8 - sourceStartIndex)
        ThrowException("Count must be > 0 and < sourceLength * 8 - sourceStartIndex (bits).");

    // Pointer to current destination byte
    UnsignedByte* curDestPtr = destination;
    // Pointer to current source byte
    UnsignedByte* sourcePtr = source;
    // Total bits to write
    size_t totalCount = count;

    // Move destination pointer to first byte
    size_t destByteIndex = destStartIndex / 8;
    curDestPtr += destByteIndex;
    // First destination bit to write to (from left to right)
    size_t destBitIndex = destStartIndex % 8;

    // Calculate direction to move in source bytes
    int step = 1;
    if (isSourceOrderReversed) {
        step = -1;
        sourcePtr = source + sourceLength - 1;
    }
    // Move source pointer to the highest byte
    sourcePtr += step * ((int)sourceStartIndex / 8);
    // Source bit index (from right to left)
    size_t sourceBitIndex = sourceStartIndex % 8;

    // Calculate number of bit to read then write for the 1st time
    size_t curCount = 8 - destBitIndex;
    size_t sourceCurCount =  8 - sourceBitIndex;
    if (curCount > sourceCurCount) {
        curCount = sourceCurCount;
    }
    if (curCount > count) {
        curCount = count;
    }

    while (totalCount > 0) {
        CopyBits(*sourcePtr, sourceBitIndex, curDestPtr, destBitIndex, curCount);
        // Should move to next source byte
        sourceBitIndex += curCount;
        if (sourceBitIndex >= 8) {
            sourceBitIndex = 0;
            sourcePtr += step;
        }
        // Should move to next destination byte
        destBitIndex += curCount;
        if (destBitIndex >= 8) {
            curDestPtr += 1;
            destBitIndex = 0;
        }
        // Calculate number of bits will be read
        totalCount -= curCount;
        curCount = totalCount > 8 ? 8 : totalCount;
        if (curCount > count) {
            curCount = count;
        }
        size_t offset = 8 - destBitIndex;
        if (curCount > offset) {
            curCount = offset;
        }
        offset = 8 - sourceBitIndex;
        if (curCount > offset) {
            curCount = offset;
        }
    }
}

/// Return array of 6 items contains QR aligment locations for version 2...40
inline const UnsignedByte * AlignmentLocations(UnsignedByte version)
{
    if (version < 2 || version > QR_MAX_VERSION) {
        ThrowException("Version must be 2...40");
    }
    static const UnsignedByte data[][6] = {
       {  18,    0,    0,    0,    0,    0},
       {  22,    0,    0,    0,    0,    0},
       {  26,    0,    0,    0,    0,    0},
       {  30,    0,    0,    0,    0,    0},
       {  34,    0,    0,    0,    0,    0},
       {  22,   38,    0,    0,    0,    0},
       {  24,   42,    0,    0,    0,    0},
       {  26,   46,    0,    0,    0,    0},
       {  28,   50,    0,    0,    0,    0},
       {  30,   54,    0,    0,    0,    0},
       {  32,   58,    0,    0,    0,    0},
       {  34,   62,    0,    0,    0,    0},
       {  26,   46,   66,    0,    0,    0},
       {  26,   48,   70,    0,    0,    0},
       {  26,   50,   74,    0,    0,    0},
       {  30,   54,   78,    0,    0,    0},
       {  30,   56,   82,    0,    0,    0},
       {  30,   58,   86,    0,    0,    0},
       {  34,   62,   90,    0,    0,    0},
       {  28,   50,   72,   94,    0,    0},
       {  26,   50,   74,   98,    0,    0},
       {  30,   54,   78,  102,    0,    0},
       {  28,   54,   80,  106,    0,    0},
       {  32,   58,   84,  110,    0,    0},
       {  30,   58,   86,  114,    0,    0},
       {  34,   62,   90,  118,    0,    0},
       {  26,   50,   74,   98,  122,    0},
       {  30,   54,   78,  102,  126,    0},
       {  26,   52,   78,  104,  130,    0},
       {  30,   56,   82,  108,  134,    0},
       {  34,   60,   86,  112,  138,    0},
       {  30,   58,   86,  114,  142,    0},
       {  34,   62,   90,  118,  146,    0},
       {  30,   54,   78,  102,  126,  150},
       {  24,   50,   76,  102,  128,  154},
       {  28,   54,   80,  106,  132,  158},
       {  32,   58,   84,  110,  136,  162},
       {  26,   54,   82,  110,  138,  166},
       {  30,   58,   86,  114,  142,  170},
    };
    return data[version - 2];
}

inline void ValidateNumeric(const UnsignedByte * data, size_t length)
{
    static const char * pattern = "0123456789";
    for (size_t index = 0; index < length; index++) {
        UnsignedByte byte = data[index];
        if (strchr(pattern, (char)byte) == NULL) {
            std::string msg = "Invalid data for Numeric mode [";
            msg.append(std::to_string(index));
            msg.append("] ");
            msg.append(std::to_string(byte));
            ThrowException(std::move(msg));
        }
    }
}

inline void ValidateAlphaNumeric(const UnsignedByte * data, size_t length)
{
    static const char * pattern = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
    for (size_t index = 0; index < length; index++) {
        UnsignedByte byte = data[index];
        if (strchr(pattern, (char)byte) == NULL) {
            std::string msg = "Invalid data for AlphaNumeric mode [";
            msg.append(std::to_string(index));
            msg.append("] ");
            msg.append(std::to_string(byte));
            ThrowException(std::move(msg));
        }
    }
}

inline void ValidateKanji(const UnsignedByte * data, size_t length) {
    // Only accept 2 bytes ShiftJIS characters
    if ((length % 2) > 0) {
        ThrowException("Invalid data for Kanji mode");
    }
    for (size_t index = 0; index < length; index += 2) {
        Unsigned2Bytes curChar = 0;
        UnsignedByte * curCharPtr = (UnsignedByte*)&curChar;
        if constexpr (common::isLittleEndian) {
            *curCharPtr = data[index + 1];
            *(curCharPtr + 1) = data[index];
        } else {
            *curCharPtr = data[index];
            *(curCharPtr + 1) = data[index + 1];
        }
        bool isValid = ((curChar >= 0x8140) && (curChar <= 0x9FFC)) ||
                       ((curChar >= 0xE040) && (curChar <= 0xEBBF));
        if (not isValid) {
            std::string msg = "Invalid data for Kanji mode [";
            msg.append(std::to_string(index));
            msg.append("] ");
            msg.append(std::to_string(data[index]));
            msg.append(" ");
            msg.append(std::to_string(data[index + 1]));
            msg.append(": ");
            msg.append(std::to_string(curChar));
            ThrowException(std::move(msg));
        }
    }
}

//todo, replace with regex
inline void ValidateInputBytes(EncodingMode mode, const UnsignedByte* data, size_t length) {
    switch (mode) {
    case EncodingMode::NUMERIC:
        ValidateNumeric(data, length);
        break;
    case EncodingMode::ALPHA_NUMERIC:
        ValidateAlphaNumeric(data, length);
        break;
    case EncodingMode::KANJI:
        ValidateKanji(data, length);
        break;
    default:
        break;
    }
}

enum class EncodingExtraMode
{
    /// none
    NONE,
    /// MicroQR
    MICRO_QR,
    /// FNC1 First Position
    FNC1_FIRST,
    /// FNC2 Second Position
    FNC1_SECOND,
};

struct Polynomial {
    size_t length{0};
    UnsignedByte* terms{nullptr};
    inline static constexpr std::array<UnsignedByte, 256> EXP = []() constexpr {
        size_t xVal = 1;
        std::array<UnsignedByte, 256> res{};
        for (size_t index = 0; index < 255; index++) {
            res[index] = static_cast<UnsignedByte>(xVal);
            if (xVal <<= 1; xVal >= 256) {
                xVal ^= 0x11D;
            }
        }
        return res;
    }();

    inline static constexpr std::array<UnsignedByte, 256> LOG = []() constexpr {
        size_t xVal = 1;
        std::array<UnsignedByte, 256> res{};
        for (size_t index = 0; index < 255; index++) {
            res[xVal] = static_cast<UnsignedByte>(index);
            if (xVal <<= 1; xVal >= 256) {
                xVal ^= 0x11D;
            }
        }
        return res;
    }();

    ~Polynomial()
    {
        if (terms) {
            delete[] terms;
            terms = nullptr;
        }
    }

    Polynomial (size_t count)
    {
        length = count;
        terms = Allocate(count);
    }

    Polynomial (Polynomial &other)
    {
        length = other.length;
        terms = new UnsignedByte [length];
        for (size_t index = 0; index < length; index += 1) {
            terms[index] = other.terms[index];
        }
    }
    
    void operator=(Polynomial other)
    {
        delete[] terms;
        length = other.length;
        terms = new UnsignedByte [length];
        for (size_t index = 0; index < length; index += 1) {
            terms[index] = other.terms[index];
        }
    }

    static UnsignedByte Multiple(UnsignedByte left, UnsignedByte right)
    {
        if (left == 0 || right == 0) {
            return 0;
        }
        return EXP[(LOG[left] + LOG[right]) % 255];
    }

    static UnsignedByte Power(UnsignedByte value, UnsignedByte power)
    {
        return EXP[(LOG[value] * power) % 255];
    }


    static Polynomial PolyMultiple(Polynomial self, Polynomial other)
    {
        Polynomial result(self.length + other.length - 1);
        for (size_t jndex = 0; jndex < other.length; jndex += 1) {
            for (size_t index = 0; index < self.length; index += 1) {
                result.terms[index + jndex] ^= Multiple(self.terms[index], other.terms[jndex]);
            }
        }
        return result;
    }

    Polynomial GetGeneratorPoly(size_t count)
    {
        Polynomial result(1);
        result.terms[0] = 1;
        for (size_t index = 0; index < count; index += 1) {
            Polynomial arg(2);
            arg.terms[0] = 1;
            arg.terms[1] = Power(2, index);
            result = PolyMultiple(result, arg);
        }
        return result;
    }

    Polynomial GetErrorCorrections(size_t count)
    {
        if (length + count > 255)
            ThrowException("Internal error: invalid message length to calculate Error Corrections");

        Polynomial gen = GetGeneratorPoly(count);
        Polynomial buffer(length + gen.length - 1);
        for (size_t index = 0; index < length; index += 1) {
            buffer.terms[index] = terms[index];
        }
        for (size_t index = 0; index < length; index += 1) {
            UnsignedByte coef = buffer.terms[index];
            if (coef != 0) {
                for (size_t jndex = 1; jndex < gen.length; jndex += 1) {
                    buffer.terms[index + jndex] ^= Multiple(gen.terms[jndex], coef);
                }
            }
        }
        Polynomial result(buffer.length - length);
        for (size_t index = length; index < buffer.length; index += 1) {
            result.terms[index - length] = buffer.terms[index];
        }
        return result;
    }
};

struct QRMatrixExtraMode//todo, remove raw pointer
{
    /// Extra mode
    EncodingExtraMode mode{EncodingExtraMode::NONE};
    /// A C-String represents Application Indicator for FNC1 Second Position.
    /// 2 valid forms:
    ///  - Single ASCII character [a-z][A-Z] eg. `a`.
    ///  - Two ditis number eg. `01`
    UnsignedByte * appIndicator{nullptr};
    UnsignedByte appIndicatorLength{0};

    ~QRMatrixExtraMode()
    {
        if (appIndicator != nullptr) {
            delete[] appIndicator;
            appIndicator = nullptr;
        }
    }

    QRMatrixExtraMode(QRMatrixExtraMode & other)
    {
        mode = other.mode;
        appIndicatorLength = other.appIndicatorLength;
        appIndicator = new UnsignedByte[appIndicatorLength];
        for (size_t index = 0; index < appIndicatorLength; index++) {
            appIndicator[index] = other.appIndicator[index];
        }
    }
    void operator=(QRMatrixExtraMode other)
    {
        if (appIndicator != nullptr) {
            delete[] appIndicator;
            appIndicator = nullptr;
        }
        mode = other.mode;
        appIndicatorLength = other.appIndicatorLength;
        appIndicator = new UnsignedByte[appIndicatorLength];
        for (size_t index = 0; index < appIndicatorLength; index++) {
            appIndicator[index] = other.appIndicator[index];
        }
    }

    /// Default init, none
    QRMatrixExtraMode() = default;
    /// Init for MicroQR or FCN1 First Position mode
    QRMatrixExtraMode(EncodingExtraMode mode) :mode(mode)
    {
    }
    /// Init for FCN2 Second Positon mode
    QRMatrixExtraMode(UnsignedByte* appId, UnsignedByte appIdLen)
    {
        mode = EncodingExtraMode::FNC1_SECOND;
        appIndicatorLength = appIdLen;
        appIndicator = new UnsignedByte[appIdLen];
        for (size_t index = 0; index < appIdLen; index++) {
            appIndicator[index] = appId[index];
        }
    }
};

struct ErrorCorrectionInfo
{
    UnsignedByte version{0};
    ErrorCorrectionLevel level{ErrorCorrectionLevel::LOW};
    Unsigned2Bytes codewords{0};
    Unsigned2Bytes ecCodewordsPerBlock{0};
    Unsigned2Bytes group1Blocks{0};
    Unsigned2Bytes group1BlockCodewords{0};
    Unsigned2Bytes group2Blocks{0};
    Unsigned2Bytes group2BlockCodewords{0};

    ErrorCorrectionInfo() = default;
    ErrorCorrectionInfo(UnsignedByte ver, ErrorCorrectionLevel level, const Unsigned2Bytes data[6])
    {
        version = ver;
        level = level;
        codewords = data[0];
        ecCodewordsPerBlock = data[1];
        group1Blocks = data[2];
        group1BlockCodewords = data[3];
        group2Blocks = data[4];
        group2BlockCodewords = data[5];
    }

    /// Number of bytes for QR data.
    /// @ref: https://www.thonky.com/qr-code-tutorial/error-correction-table.
    static ErrorCorrectionInfo GetErrorCorrectionInfo (UnsignedByte version, ErrorCorrectionLevel level)
    {
        static const Unsigned2Bytes map[40][4][6] = {
            {{  19,  7,  1,  19,  0,   0}, {  16, 10,  1,  16,  0,   0}, {  13, 13,  1,  13,  0,   0}, {   9, 17,  1,   9,  0,   0}},
            {{  34, 10,  1,  34,  0,   0}, {  28, 16,  1,  28,  0,   0}, {  22, 22,  1,  22,  0,   0}, {  16, 28,  1,  16,  0,   0}},
            {{  55, 15,  1,  55,  0,   0}, {  44, 26,  1,  44,  0,   0}, {  34, 18,  2,  17,  0,   0}, {  26, 22,  2,  13,  0,   0}},
            {{  80, 20,  1,  80,  0,   0}, {  64, 18,  2,  32,  0,   0}, {  48, 26,  2,  24,  0,   0}, {  36, 16,  4,   9,  0,   0}},
            {{ 108, 26,  1, 108,  0,   0}, {  86, 24,  2,  43,  0,   0}, {  62, 18,  2,  15,  2,  16}, {  46, 22,  2,  11,  2,  12}},
            {{ 136, 18,  2,  68,  0,   0}, { 108, 16,  4,  27,  0,   0}, {  76, 24,  4,  19,  0,   0}, {  60, 28,  4,  15,  0,   0}},
            {{ 156, 20,  2,  78,  0,   0}, { 124, 18,  4,  31,  0,   0}, {  88, 18,  2,  14,  4,  15}, {  66, 26,  4,  13,  1,  14}},
            {{ 194, 24,  2,  97,  0,   0}, { 154, 22,  2,  38,  2,  39}, { 110, 22,  4,  18,  2,  19}, {  86, 26,  4,  14,  2,  15}},
            {{ 232, 30,  2, 116,  0,   0}, { 182, 22,  3,  36,  2,  37}, { 132, 20,  4,  16,  4,  17}, { 100, 24,  4,  12,  4,  13}},
            {{ 274, 18,  2,  68,  2,  69}, { 216, 26,  4,  43,  1,  44}, { 154, 24,  6,  19,  2,  20}, { 122, 28,  6,  15,  2,  16}},
            {{ 324, 20,  4,  81,  0,   0}, { 254, 30,  1,  50,  4,  51}, { 180, 28,  4,  22,  4,  23}, { 140, 24,  3,  12,  8,  13}},
            {{ 370, 24,  2,  92,  2,  93}, { 290, 22,  6,  36,  2,  37}, { 206, 26,  4,  20,  6,  21}, { 158, 28,  7,  14,  4,  15}},
            {{ 428, 26,  4, 107,  0,   0}, { 334, 22,  8,  37,  1,  38}, { 244, 24,  8,  20,  4,  21}, { 180, 22, 12,  11,  4,  12}},
            {{ 461, 30,  3, 115,  1, 116}, { 365, 24,  4,  40,  5,  41}, { 261, 20, 11,  16,  5,  17}, { 197, 24, 11,  12,  5,  13}},
            {{ 523, 22,  5,  87,  1,  88}, { 415, 24,  5,  41,  5,  42}, { 295, 30,  5,  24,  7,  25}, { 223, 24, 11,  12,  7,  13}},
            {{ 589, 24,  5,  98,  1,  99}, { 453, 28,  7,  45,  3,  46}, { 325, 24, 15,  19,  2,  20}, { 253, 30,  3,  15, 13,  16}},
            {{ 647, 28,  1, 107,  5, 108}, { 507, 28, 10,  46,  1,  47}, { 367, 28,  1,  22, 15,  23}, { 283, 28,  2,  14, 17,  15}},
            {{ 721, 30,  5, 120,  1, 121}, { 563, 26,  9,  43,  4,  44}, { 397, 28, 17,  22,  1,  23}, { 313, 28,  2,  14, 19,  15}},
            {{ 795, 28,  3, 113,  4, 114}, { 627, 26,  3,  44, 11,  45}, { 445, 26, 17,  21,  4,  22}, { 341, 26,  9,  13, 16,  14}},
            {{ 861, 28,  3, 107,  5, 108}, { 669, 26,  3,  41, 13,  42}, { 485, 30, 15,  24,  5,  25}, { 385, 28, 15,  15, 10,  16}},
            {{ 932, 28,  4, 116,  4, 117}, { 714, 26, 17,  42,  0,   0}, { 512, 28, 17,  22,  6,  23}, { 406, 30, 19,  16,  6,  17}},
            {{1006, 28,  2, 111,  7, 112}, { 782, 28, 17,  46,  0,   0}, { 568, 30,  7,  24, 16,  25}, { 442, 24, 34,  13,  0,   0}},
            {{1094, 30,  4, 121,  5, 122}, { 860, 28,  4,  47, 14,  48}, { 614, 30, 11,  24, 14,  25}, { 464, 30, 16,  15, 14,  16}},
            {{1174, 30,  6, 117,  4, 118}, { 914, 28,  6,  45, 14,  46}, { 664, 30, 11,  24, 16,  25}, { 514, 30, 30,  16,  2,  17}},
            {{1276, 26,  8, 106,  4, 107}, {1000, 28,  8,  47, 13,  48}, { 718, 30,  7,  24, 22,  25}, { 538, 30, 22,  15, 13,  16}},
            {{1370, 28, 10, 114,  2, 115}, {1062, 28, 19,  46,  4,  47}, { 754, 28, 28,  22,  6,  23}, { 596, 30, 33,  16,  4,  17}},
            {{1468, 30,  8, 122,  4, 123}, {1128, 28, 22,  45,  3,  46}, { 808, 30,  8,  23, 26,  24}, { 628, 30, 12,  15, 28,  16}},
            {{1531, 30,  3, 117, 10, 118}, {1193, 28,  3,  45, 23,  46}, { 871, 30,  4,  24, 31,  25}, { 661, 30, 11,  15, 31,  16}},
            {{1631, 30,  7, 116,  7, 117}, {1267, 28, 21,  45,  7,  46}, { 911, 30,  1,  23, 37,  24}, { 701, 30, 19,  15, 26,  16}},
            {{1735, 30,  5, 115, 10, 116}, {1373, 28, 19,  47, 10,  48}, { 985, 30, 15,  24, 25,  25}, { 745, 30, 23,  15, 25,  16}},
            {{1843, 30, 13, 115,  3, 116}, {1455, 28,  2,  46, 29,  47}, {1033, 30, 42,  24,  1,  25}, { 793, 30, 23,  15, 28,  16}},
            {{1955, 30, 17, 115,  0,   0}, {1541, 28, 10,  46, 23,  47}, {1115, 30, 10,  24, 35,  25}, { 845, 30, 19,  15, 35,  16}},
            {{2071, 30, 17, 115,  1, 116}, {1631, 28, 14,  46, 21,  47}, {1171, 30, 29,  24, 19,  25}, { 901, 30, 11,  15, 46,  16}},
            {{2191, 30, 13, 115,  6, 116}, {1725, 28, 14,  46, 23,  47}, {1231, 30, 44,  24,  7,  25}, { 961, 30, 59,  16,  1,  17}},
            {{2306, 30, 12, 121,  7, 122}, {1812, 28, 12,  47, 26,  48}, {1286, 30, 39,  24, 14,  25}, { 986, 30, 22,  15, 41,  16}},
            {{2434, 30,  6, 121, 14, 122}, {1914, 28,  6,  47, 34,  48}, {1354, 30, 46,  24, 10,  25}, {1054, 30,  2,  15, 64,  16}},
            {{2566, 30, 17, 122,  4, 123}, {1992, 28, 29,  46, 14,  47}, {1426, 30, 49,  24, 10,  25}, {1096, 30, 24,  15, 46,  16}},
            {{2702, 30,  4, 122, 18, 123}, {2102, 28, 13,  46, 32,  47}, {1502, 30, 48,  24, 14,  25}, {1142, 30, 42,  15, 32,  16}},
            {{2812, 30, 20, 117,  4, 118}, {2216, 28, 40,  47,  7,  48}, {1582, 30, 43,  24, 22,  25}, {1222, 30, 10,  15, 67,  16}},
            {{2956, 30, 19, 118,  6, 119}, {2334, 28, 18,  47, 31,  48}, {1666, 30, 34,  24, 34,  25}, {1276, 30, 20,  15, 61,  16}},
        };
        int versionIndex = version - 1;
        if (versionIndex < 0 || versionIndex >= QR_MAX_VERSION)
            return ErrorCorrectionInfo();
        int levelIndex = ConvertErrorCorrectionLevelToIndex(level);
        const Unsigned2Bytes * data = map[versionIndex][levelIndex];
        return ErrorCorrectionInfo(version, level, data);
    }

    static ErrorCorrectionInfo GetMicroErrorCorrectionInfo(UnsignedByte version, ErrorCorrectionLevel level)
    {
        static const Unsigned2Bytes map[4][4][2] = {
            {{  3, 2}, {  3,  2}, {  3,  2}, { 3, 2}},
            {{  5, 5}, {  4,  6}, {  0,  0}, { 0, 0}},
            {{ 11, 6}, {  9,  8}, {  0,  0}, { 0, 0}},
            {{ 16, 8}, { 14, 10}, { 10, 13}, { 0, 0}},
        };
        int versionIndex = version - 1;
        if (versionIndex < 0 || versionIndex >= MICROQR_MAX_VERSION)
            return ErrorCorrectionInfo();
        int levelIndex = ConvertErrorCorrectionLevelToIndex(level);
        const Unsigned2Bytes* data = map[versionIndex][levelIndex];
        Unsigned2Bytes buffer[6];
        buffer[0] = data[0];
        buffer[1] = data[1];
        buffer[2] = 1;
        buffer[3] = data[0];
        buffer[4] = 0;
        buffer[5] = 0;
        return ErrorCorrectionInfo(version, level, buffer);
    }

    UnsignedByte GroupCount() { return group2Blocks == 0 ? 1 : 2; }
    UnsignedByte EcBlockTotalCount() { return group1Blocks + group2Blocks; }
    size_t EcCodewordsTotalCount() { return EcBlockTotalCount() * ecCodewordsPerBlock; }
};


 /// Hold information of input data for 1 QR segment.
class QRMatrixSegment // todo, refine the raw pointer
{
public:
    /// Create QR segment
    QRMatrixSegment(
        /// Encoding mode
        EncodingMode mode,
        /// Bytes sequence to encode.
        /// If mode is `Numeric`, `data` must contain only ASCII bytes of characters `[0...9]`
        /// (or other text encoding if compatiple eg. UTF-8).
        /// If mode is `AlphaNumeric, `data` must contain only ASCII bytes of characters `[0...9][A...Z] $%*+-./:`
        /// (or other text encoding if compatiple eg. UTF -8).
        /// If mode is `Kanji`, `data` must contain only 2-bytes ShiftJIS characters
        /// (each 2-bytes must be in range [0x8140...0x9FFC] & [0xE040...0xEBBF]).
        const UnsignedByte* data,
        /// Number of `data` bytes
        size_t length,
        /// Enable ECI mode with given ECI Indicator (ECI Assigment value)
        size_t eciIndicator = DEFAULT_ECI_ASSIGMENT_VALUE
    ) : m_mode(mode), m_length(length), m_eci(eciIndicator)
    {
        ValidateInputBytes(mode, data, m_length);
        if (m_length > 0) {
            m_data = Allocate(m_length);
            for (size_t index = 0; index < m_length; index++) {
                m_data[index] = data[index];
            }
        }
    }
    /// Create empty and set later
    QRMatrixSegment() = default;
    QRMatrixSegment(QRMatrixSegment &other)
    {
        m_length = other.m_length;
        m_mode = other.m_mode;
        m_eci = other.m_eci;
        if (m_length > 0) {
            m_data = Allocate(m_length);
            for (size_t index = 0; index < m_length; index++)
                m_data[index] = other.m_data[index];
        }
    }

    ~QRMatrixSegment()
    {
        if (m_data) {
            delete[] m_data;
            m_data = nullptr;
        }
    }
    /// Fill segment with given data
    void Fill(EncodingMode mode, const UnsignedByte* data, size_t length, size_t eciIndicator = DEFAULT_ECI_ASSIGMENT_VALUE)
    {
        ValidateInputBytes(mode, data, length);
        if (m_data) {
            delete[] m_data;
            m_data = nullptr;
        }
        m_mode = mode;
        m_length = length;
        m_eci = eciIndicator;
         if (m_length > 0) {
            m_data = Allocate(m_length);
            for (size_t index = 0; index < m_length; index++)
                m_data[index] = data[index];
        }
        else {
            m_data = nullptr;
        }
    }

    void operator=(QRMatrixSegment other) //todo, remove
    {
        if (m_data) {
            delete[] m_data;
            m_data = nullptr;
        }
        m_length = other.m_length;
        m_mode = other.m_mode;
        m_eci = other.m_eci;
         if (m_length > 0) {
            m_data = Allocate(m_length);
            for (size_t index = 0; index < m_length; index++)
                m_data[index] = other.m_data[index];
        }
        else {
            m_data = nullptr;
        }
    }

    EncodingMode Mode() { return m_mode; }
    size_t Length() { return m_length; }
    size_t ECI() { return m_eci; }
    UnsignedByte * Data() { return m_data; }
    /// 0 to ignore ECI Indicator.
    /// Default QR ECI indicator is 3, so we ignore too.
    /// MicroQR does not have ECI mode, so we ignore this in MicroQR.
    bool isEciHeaderRequired() { return m_eci != DEFAULT_ECI_ASSIGMENT_VALUE; }

private:
    UnsignedByte* m_data{nullptr};
    EncodingMode m_mode{EncodingMode::BYTE};
    size_t m_length{0};
    size_t m_eci{DEFAULT_ECI_ASSIGMENT_VALUE};
};

class QRMatrixBoard 
{
public:
    QRMatrixBoard() = default;
    ~QRMatrixBoard()
    {
        if (m_dimension != 0) {
            for (UnsignedByte index = 0; index < m_dimension; index += 1) {
                delete[] m_buffer[index];
            }
            delete[] m_buffer;
        }
    }

    QRMatrixBoard(QRMatrixBoard & other)
    {
        m_dimension = other.m_dimension;
        m_buffer = new UnsignedByte* [m_dimension];
        for (size_t index = 0; index < m_dimension; index ++) {
            m_buffer[index] = new UnsignedByte [m_dimension];
            for (size_t jndex = 0; jndex < m_dimension; jndex ++) {
                m_buffer[index][jndex] = other.m_buffer[index][jndex];
            }
        }
    }
    void operator=(QRMatrixBoard other)
    {
        if (m_dimension > 0) {
            for (size_t index = 0; index < m_dimension; index ++) {
                delete[] m_buffer[index];
            }
            delete[] m_buffer;
            m_buffer = nullptr;
        }
        m_dimension = other.m_dimension;
        if (m_dimension > 0) {
            m_buffer = new UnsignedByte* [m_dimension];
            for (size_t index = 0; index < m_dimension; index ++) {
                m_buffer[index] = new UnsignedByte [m_dimension];
                for (size_t jndex = 0; jndex < m_dimension; jndex ++) {
                    m_buffer[index][jndex] = other.m_buffer[index][jndex];
                }
            }
        }
    }

    /// Place holder. Internal purpose. Do not use.
    /// To create QR board, refer `QRMatrixEncoder`.
    /// This constructor is for internal purpose.
    QRMatrixBoard(UnsignedByte* data, UnsignedByte* errorCorrection, ErrorCorrectionInfo ecInfo, UnsignedByte maskId, bool isMicro)
    {
        m_dimension = isMicro ? GetMicroDimensionByVersion(ecInfo.version) : GetDimensionByVersion(ecInfo.version);
        m_buffer = new UnsignedByte* [m_dimension];
        for (UnsignedByte index = 0; index < m_dimension; index += 1) {
            m_buffer[index] = new UnsignedByte [m_dimension];
            for (UnsignedByte jndex = 0; jndex < m_dimension; jndex += 1) {
                m_buffer[index][jndex] = BoardCell::NEUTRAL;
            }
        }
        AddFinderPatterns(isMicro);
        AddSeparators(isMicro);
        if (not isMicro)
            AddAlignmentPatterns(ecInfo);
        AddTimingPatterns(isMicro);
        if (isMicro)
            AddMicroReservedAreas(ecInfo);
        else
            AddDarkAndReservedAreas(ecInfo);

        PlaceData(data, errorCorrection, ecInfo,
                    isMicro ? 0 : RemainderBitsLength(ecInfo.version), isMicro);
        UnsignedByte lastMaskId = Evaluate(maskId, isMicro);
        if (isMicro) {
            PlaceMicroFormat(lastMaskId, ecInfo);
        } else {
            PlaceFormatAndVersion(lastMaskId, ecInfo);
        }
    }

    /// Size (dimension - number of cells on each side)
    inline UnsignedByte Dimension() const { return m_dimension; }
    /// UnsignedByte[row][column]: 1 byte per cell.
    /// High word (left 4 bits) is cell type.
    /// Low word (right 4 bits) is black or white.
    /// See `BoardCell` for values.
    inline UnsignedByte ** Buffer() const { return m_buffer; }
    
    /// To print to Console
    std::string Description(bool isTypeVisible = false)
    {
        if (m_dimension == 0) return std::string{};

        std::string result = "  ";
        for (UnsignedByte index = 0; index < m_dimension; index += 1) {
            if (index % 2 > 0) {
                result.append("..");
                continue;
            }
            if (index < 10) {
                result.append("0");
            }
            result.append(std::to_string(index));
        }
        result.append("\n");
        for (UnsignedByte index = 0; index < m_dimension; index += 1) {
            if (index < 10) {
                result.append("0");
            }
            result.append(std::to_string(index));
            for (UnsignedByte jndex = 0; jndex < m_dimension; jndex += 1) {
                UnsignedByte byte = m_buffer[index][jndex];
                if (!isTypeVisible) {
                    UnsignedByte low = byte & BoardCell::LOWMASK;
                    if (low == BoardCell::UNSET) {
                        result.append(" □");
                    } else if (low == BoardCell::SET) {
                        result.append(" ■");
                    } else {
                        result.append("  ");
                    }
                } else {
                    if (byte == BoardCell::NEUTRAL) {
                        result.append(" +");
                    } else {
                        UnsignedByte low = byte & BoardCell::LOWMASK;
                        UnsignedByte high = byte & BoardCell::HIGHMASK;
                        bool isFunc = (byte & BoardCell::FUNCMASK) > 0;
                        if (high == BoardCell::FORMAT || high == BoardCell::VERSION) {
                            result.append("®");
                        } else if (isFunc) {
                            result.append("•");
                        } else if (high == BoardCell::ERROR_CORRECTION) {
                            result.append(".");
                        } else {
                            result.append(" ");
                        }
                        switch (low) {
                        case BoardCell::UNSET:
                            result.append("□");
                            break;
                        case BoardCell::SET:
                            result.append("■");
                            break;
                        default:
                            if (high == BoardCell::FORMAT || high == BoardCell::VERSION) {
                                result.append("®");
                            } else {
                                result.append(" ");
                            }
                            break;
                        }
                    }
                }
            }
            result.append("\n");
        }
        return result;
    }
private:
    void AddFinderPatterns(bool isMicro) 
    {
        AddFinderPattern(0, 0);
        if (isMicro) return;
        AddFinderPattern(0, m_dimension - 7);
        AddFinderPattern(m_dimension - 7, 0);
    }
    
    void AddFinderPattern(UnsignedByte row, UnsignedByte column)
    {
        SetSquare(row, column, 7, false, true, BoardCell::FINDER);
        SetSquare(row + 1, column + 1, 5, false, false, BoardCell::FINDER);
        SetSquare(row + 2, column + 2, 3, true, true, BoardCell::FINDER);
    }

    void SetSquare(UnsignedByte row, UnsignedByte column, UnsignedByte size, bool isFill, bool isSet, UnsignedByte prefix)
    {
        UnsignedByte value = (isSet ? BoardCell::SET : BoardCell::UNSET) | prefix;
        if (isFill) {
            for (UnsignedByte rIndex = row; rIndex < row + size; rIndex += 1) {
                for (UnsignedByte cIndex = column; cIndex < column + size; cIndex += 1) {
                    m_buffer[rIndex][cIndex] = value;
                }
            }
        } else {
            for (UnsignedByte index = 0; index < size; index += 1) {
                m_buffer[row + index][column] = value;
                m_buffer[row][column + index] = value;
                m_buffer[row + index][column + size - 1] = value;
                m_buffer[row + size - 1][column + index] = value;
            }
        }
    }

    void AddSeparators(bool isMicro)
    {
        UnsignedByte value = BoardCell::UNSET | BoardCell::SEPARATOR;
        if (isMicro) {
            for (UnsignedByte index = 0; index < 8; index += 1) {
                m_buffer[7][index] = value;
                m_buffer[index][7] = value;
            }
        } else {
            for (UnsignedByte index = 0; index < 8; index += 1) {
                m_buffer[7][index] = value;
                m_buffer[index][7] = value;

                m_buffer[7][m_dimension - index - 1] = value;
                m_buffer[index][m_dimension - 8] = value;

                m_buffer[m_dimension - 8][index] = value;
                m_buffer[m_dimension - index - 1][7] = value;
            }
        }
    }

    void AddAlignmentPatterns(ErrorCorrectionInfo ecInfo)
    {
        if (ecInfo.version < 2) return;
        const UnsignedByte* array = AlignmentLocations(ecInfo.version);
        AddAlignmentPattern(6, 6);
        for (UnsignedByte index = 0; index < 6; index += 1) {
            UnsignedByte value = array[index];
            if (value > 0) {
                AddAlignmentPattern(6, value);
                AddAlignmentPattern(value, 6);
                AddAlignmentPattern(value, value);
                for (UnsignedByte jndex = 0; jndex < 6; jndex += 1) {
                    UnsignedByte value2 = array[jndex];
                    if (value2 > 0 && value != value2) {
                        AddAlignmentPattern(value, value2);
                        AddAlignmentPattern(value2, value);
                    }
                }
            }
        }
    }

    void AddAlignmentPattern(UnsignedByte row, UnsignedByte column)
    {
        UnsignedByte tlRow = row - 2;
        UnsignedByte tlCol = column - 2;
        for (UnsignedByte rIndex = tlRow; rIndex < tlRow + 5; rIndex += 1) {
            for (UnsignedByte cIndex = tlCol; cIndex < tlCol + 5; cIndex += 1) {
                if (m_buffer[rIndex][cIndex] != BoardCell::NEUTRAL) return;
            }
        }
        SetSquare(tlRow, tlCol, 5, false, true, BoardCell::ALIGNMENT);
        SetSquare(tlRow + 1, tlCol + 1, 3, false, false, BoardCell::ALIGNMENT);
        m_buffer[row][column] = BoardCell::SET | BoardCell::ALIGNMENT;
    }

    void AddTimingPatterns(bool isMicro)
    {
        UnsignedByte valueSet = BoardCell::TIMING | BoardCell::SET;
        UnsignedByte valueUnset = BoardCell::TIMING | BoardCell::UNSET;
        UnsignedByte offset = isMicro ? 0 : 6;
        for (UnsignedByte index = 6; index < m_dimension - offset; index += 1) {
            if ((index % 2) == 0) {
                m_buffer[offset][index] = valueSet;
                m_buffer[index][offset] = valueSet;
            } else {
                m_buffer[offset][index] = valueUnset;
                m_buffer[index][offset] = valueUnset;
            }
        }
    }

    void AddMicroReservedAreas(ErrorCorrectionInfo ecInfo)
    {
        for (UnsignedByte index = 0; index < 8; index += 1) {
            if (m_buffer[8][index] == BoardCell::NEUTRAL) {
                m_buffer[8][index] = BoardCell::FORMAT | BoardCell::UNSET;
            }
            if (m_buffer[index][8] == BoardCell::NEUTRAL) {
                m_buffer[index][8] = BoardCell::FORMAT | BoardCell::UNSET;
            }
        }
        m_buffer[8][8] = BoardCell::FORMAT | BoardCell::UNSET;
    }

    void AddDarkAndReservedAreas(ErrorCorrectionInfo ecInfo)
    {
        // Dark cell
        m_buffer[m_dimension - 8][8] = BoardCell::DARK | BoardCell::SET;
        // Reseved cells for format
        for (UnsignedByte index = 0; index < 8; index += 1) {
            if (m_buffer[8][index] == BoardCell::NEUTRAL) {
                m_buffer[8][index] = BoardCell::FORMAT | BoardCell::UNSET;
            }
            if (m_buffer[index][8] == BoardCell::NEUTRAL) {
                m_buffer[index][8] = BoardCell::FORMAT | BoardCell::UNSET;
            }
            if (m_buffer[m_dimension - index - 1][8] == BoardCell::NEUTRAL) {
                m_buffer[m_dimension - index - 1][8] = BoardCell::FORMAT | BoardCell::UNSET;
            }
            m_buffer[8][m_dimension - index - 1] = BoardCell::FORMAT | BoardCell::UNSET;
        }
        m_buffer[8][8] = BoardCell::FORMAT | BoardCell::UNSET;
        if (ecInfo.version < 7) return;
        for (UnsignedByte index = 0; index < 3; index += 1) {
            for (UnsignedByte jndex = 0; jndex < 6; jndex += 1) {
                m_buffer[jndex][m_dimension - 9 - index] = BoardCell::VERSION | BoardCell::UNSET;
                m_buffer[m_dimension - 9 - index][jndex] = BoardCell::VERSION | BoardCell::UNSET;
            }
        }
    }

    /// Fill content, EC data into QR board
    void PlaceData(UnsignedByte* data, UnsignedByte* errorCorrection, ErrorCorrectionInfo ecInfo, UnsignedByte remainderCount, bool isMicro)
    {
        bool isMicroV13 = isMicro && (ecInfo.version == 1 || ecInfo.version == 3);
        bool isUpward = true;
        int column = m_dimension - 1;
        size_t byteIndex = 0;
        UnsignedByte bitIndex = 0;
        size_t bitCount = 0;
        size_t dataBitTotal = ecInfo.codewords * 8;
        if (isMicroV13) {
            dataBitTotal -= 4;
        }
        size_t ecBitTotal = ecInfo.EcCodewordsTotalCount() * 8;
        UnsignedByte phase = 0;
        bool isCompleted = false;

        while (column >= 0 && !isCompleted) {
            int startValue = 0;
            int endValue = m_dimension - 1;
            int step = 1;
            if (isUpward) {
                startValue = m_dimension - 1;
                endValue = 0;
                step = -1;
            }
            for (int row = startValue; (isUpward ? row >= endValue : row <= endValue) && !isCompleted; row += step) {
                if (m_buffer[row][column] == BoardCell::NEUTRAL) {
                    FillDataBit(data, errorCorrection, phase,
                                &byteIndex, &bitIndex, row, column);
                    isCompleted = CheckFilledBit(&bitCount, dataBitTotal, ecBitTotal, &byteIndex,
                                                 &bitIndex, &phase, remainderCount);
                }
                if (column > 0 && m_buffer[row][column - 1] == BoardCell::NEUTRAL) {
                    FillDataBit(data, errorCorrection, phase,
                                &byteIndex, &bitIndex, row, column - 1);
                    isCompleted = CheckFilledBit(
                        &bitCount, dataBitTotal, ecBitTotal, &byteIndex,
                        &bitIndex, &phase, remainderCount
                        );
                }
            }
            column -= 2;
            if (!isMicro && column == 6) {
                column -= 1;
            }
            isUpward = !isUpward;
        }
    }

    /// Fill data bit into cell; increase data byte & bit index.
    void FillDataBit(UnsignedByte* data, UnsignedByte* errorCorrection, UnsignedByte phase, size_t* byteIndex, UnsignedByte* bitIndex, int row, int column)
    {
        UnsignedByte value = 0;
        UnsignedByte mask = 0b10000000 >> *bitIndex;
        UnsignedByte prefix = 0x00;

        switch (phase) {
        case 0:
            value = data[*byteIndex];
            break;
        case 1:
            value = errorCorrection[*byteIndex];
            prefix = BoardCell::ERROR_CORRECTION;
            break;
        case 2:
            value = 0;
            mask = 0;
            prefix = BoardCell::REMAINDER;
            break;
        default:
            break;
        }
        if ((value & mask) > 0) {
            m_buffer[row][column] = BoardCell::SET | prefix;
        } else {
            m_buffer[row][column] = BoardCell::UNSET | prefix;
        }
        *bitIndex += 1;
        if (*bitIndex > 7) {
            *bitIndex = 0;
            *byteIndex += 1;
        }
    }

    bool CheckFilledBit(size_t* bitCount, size_t dataBitTotal, size_t ecBitTotal, size_t* byteIndex, UnsignedByte* bitIndex, UnsignedByte* phase, UnsignedByte remainderCount)
    {
        *bitCount += 1;
        switch (*phase) {
        case 0:
            if (*bitCount >= dataBitTotal) {
                *bitCount = 0;
                *phase += 1;
                *byteIndex = 0;
                *bitIndex = 0;
            }
            break;
        case 1:
            if (*bitCount >= ecBitTotal) {
                *bitCount = 0;
                *phase += 1;
                *byteIndex = 0;
                *bitIndex = 0;
                if (remainderCount == 0) {
                    return true;
                }
            }
            break;
        case 2:
            return *bitCount >= remainderCount;
        default:
            break;
        }
        return false;
    }

    static UnsignedByte RemainderBitsLength(UnsignedByte version)
    {
        if (version >= 2 && version <= 6) {
            return 7;
        }
        if ((version >= 14 && version <= 20) || (version >= 28 && version <= 34)) {
            return 3;
        }
        if (version >= 21 && version <= 27) {
            return 4;
        }
        return 0;
    }

    UnsignedByte** Mask(UnsignedByte maskNum) 
    {
        UnsignedByte** result = new UnsignedByte* [m_dimension];
        for (UnsignedByte row = 0; row < m_dimension; row += 1) {
            result[row] = new UnsignedByte [m_dimension];
            for (UnsignedByte column = 0; column < m_dimension; column += 1) {
                UnsignedByte byte = m_buffer[row][column];
                UnsignedByte low = byte & BoardCell::LOWMASK;
                UnsignedByte high = byte & BoardCell::HIGHMASK;
                bool isFunc = (byte & BoardCell::FUNCMASK) > 0;
                if (isFunc) {
                    result[row][column] = byte;
                    continue;
                }
                if (low == BoardCell::SET) {
                    low = BoardCell::UNSET;
                } else if (low == BoardCell::UNSET) {
                    low = BoardCell::SET;
                }
                UnsignedByte maskedByte = low | high;
                switch (maskNum) {
                case 0:
                    if (((row + column) % 2) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                case 1:
                    if ((row % 2) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                case 2:
                    if ((column % 3) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                case 3:
                    if (((row + column) % 3) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                case 4:
                    if (((UnsignedByte)(floor((double)row / 2) + floor((double)column / 3)) % 2) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                case 5:
                    if (((row * column) % 2 + (row * column) % 3) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                case 6:
                    if ((((row * column) % 2 + (row * column) % 3) % 2) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                case 7:
                    if ((((row + column) % 2 + (row * column) % 3) % 2) == 0) {
                        result[row][column] = maskedByte;
                    } else {
                        result[row][column] = byte;
                    }
                    break;
                }
            }
        }
        return result;
    }

    size_t EvaluateCondition1(UnsignedByte** maskedBoard)
    {
        size_t result = 0;
        UnsignedByte sameColorCount = 0;
        UnsignedByte curColor = BoardCell::DARK;
        for (UnsignedByte row = 0; row < m_dimension; row += 1) {
            for (UnsignedByte column = 0; column < m_dimension; column += 1) {
                UnsignedByte cell = maskedBoard[row][column] & BoardCell::LOWMASK;
                if (cell == curColor) {
                    sameColorCount += 1;
                } else {
                    curColor = cell;
                    if (sameColorCount >= 5) {
                        result += 3 + (sameColorCount - 5);
                    }
                    sameColorCount = 0;
                }
            }
        }
        if (sameColorCount >= 5) {
            result += 3 + (sameColorCount - 5);
        }
        sameColorCount = 0;
        for (UnsignedByte column = 0; column < m_dimension; column += 1) {
            for (UnsignedByte row = 0; row < m_dimension; row += 1) {
                UnsignedByte cell = maskedBoard[row][column] & BoardCell::LOWMASK;
                if (cell == curColor) {
                    sameColorCount += 1;
                } else {
                    curColor = cell;
                    if (sameColorCount >= 5) {
                        result += 3 + (sameColorCount - 5);
                    }
                    sameColorCount = 0;
                }
            }
        }
        if (sameColorCount >= 5) {
            result += 3 + (sameColorCount - 5);
        }
        return result;
    }

    size_t EvaluateCondition2(UnsignedByte** maskedBoard)
    {
        size_t result = 0;
        for (UnsignedByte row = 0; row < m_dimension - 1; row += 1) {
            for (UnsignedByte column = 0; column < m_dimension - 1; column += 1) {
                UnsignedByte cell = maskedBoard[row][column] & BoardCell::LOWMASK;
                UnsignedByte cell1 = maskedBoard[row][column + 1] & BoardCell::LOWMASK;
                UnsignedByte cell2 = maskedBoard[row + 1][column] & BoardCell::LOWMASK;
                UnsignedByte cell3 = maskedBoard[row + 1][column + 1] & BoardCell::LOWMASK;
                if (cell == cell1 && cell == cell2 && cell == cell3) {
                    result += 3;
                }
            }
        }
        return result;
    }

    size_t EvaluateCondition3(UnsignedByte** maskedBoard)
    {
        static const UnsignedByte pattern1[11] = {
            BoardCell::SET,
            BoardCell::UNSET,
            BoardCell::SET,
            BoardCell::SET,
            BoardCell::SET,
            BoardCell::UNSET,
            BoardCell::SET,
            BoardCell::UNSET,
            BoardCell::UNSET,
            BoardCell::UNSET,
            BoardCell::UNSET
        };
        static const UnsignedByte pattern2[11] = {
            BoardCell::UNSET,
            BoardCell::UNSET,
            BoardCell::UNSET,
            BoardCell::UNSET,
            BoardCell::SET,
            BoardCell::UNSET,
            BoardCell::SET,
            BoardCell::SET,
            BoardCell::SET,
            BoardCell::UNSET,
            BoardCell::SET
        };
        size_t result = 0;
        for (UnsignedByte row = 0; row < m_dimension; row += 1) {
            for (UnsignedByte column = 0; column < m_dimension - 11; column += 1) {
                bool isMatched1 = true;
                bool isMatched2 = true;
                for (UnsignedByte index = 0; index < 11; index += 1) {
                    UnsignedByte cell = maskedBoard[row][column + index] & BoardCell::LOWMASK;
                    if (cell != pattern1[index]) {
                        isMatched1 = false;
                    }
                    if (cell != pattern2[index]) {
                        isMatched2 = false;
                    }
                    if (!isMatched1 && !isMatched2) {
                        break;
                    }
                }
                if (isMatched1 || isMatched2) {
                    result += 40;
                }
            }
        }
        for (UnsignedByte column = 0; column < m_dimension; column += 1) {
            for (UnsignedByte row = 0; row < m_dimension - 11; row += 1) {
                bool isMatched1 = true;
                bool isMatched2 = true;
                for (UnsignedByte index = 0; index < 11; index += 1) {
                    UnsignedByte cell = maskedBoard[row + index][column] & BoardCell::LOWMASK;
                    if (cell != pattern1[index]) {
                        isMatched1 = false;
                    }
                    if (cell != pattern2[index]) {
                        isMatched2 = false;
                    }
                    if (!isMatched1 && !isMatched2) {
                        break;
                    }
                }
                if (isMatched1 || isMatched2) {
                    result += 40;
                }
            }
        }
        return result;
    }

    size_t EvaluateCondition4(UnsignedByte** maskedBoard) {
        size_t total = m_dimension * m_dimension;
        size_t darkCount = 0;
        for (UnsignedByte row = 0; row < m_dimension; row += 1) {
            for (UnsignedByte column = 0; column < m_dimension; column += 1) {
                UnsignedByte cell = maskedBoard[row][column] & BoardCell::LOWMASK;
                if (cell == BoardCell::SET) {
                    darkCount += 1;
                }
            }
        }
        double percent = ((double)darkCount / (double)total) * 100.0f;
        size_t pre5 = percent / 5;
        pre5 *= 5;
        size_t next5 = pre5 + 5;
        pre5 = abs((int)pre5 - 50);
        next5 = abs((int)next5 - 50);
        pre5 = pre5 / 5;
        next5 = next5 / 5;

        size_t result = std::min(pre5, next5) * 10;
        return result;
    }

    size_t EvaluateMicro(UnsignedByte** maskedBoard)
    {
        UnsignedByte sum1 = 0;
        UnsignedByte sum2 = 0;
        for (size_t index = 0; index < m_dimension; index += 1) {
            UnsignedByte cell = maskedBoard[index][m_dimension - 1] & BoardCell::LOWMASK;
            if (cell == BoardCell::SET) {
                sum1 += 1;
            }
            cell = maskedBoard[m_dimension - 1][index] & BoardCell::LOWMASK;
            if (cell == BoardCell::SET) {
                sum2 += 1;
            }
        }
        if (sum1 <= sum2) {
            return sum1 * 16 + sum2;
        }
        return sum2 * 16 + sum1;
    }

    UnsignedByte Evaluate(UnsignedByte maskId, bool isMicro)
    {
        static UnsignedByte microMaskIdMap[4] = {1, 4, 6, 7};
        bool isCustomMask = isMicro ? (maskId < 4) : (maskId < 8);
        if (isCustomMask) {
            UnsignedByte mId = isMicro ? microMaskIdMap[maskId] : maskId;
            UnsignedByte** mBoard = Mask(mId);
            for (UnsignedByte row = 0; row < m_dimension; row += 1) {
                for (UnsignedByte column = 0; column < m_dimension; column += 1) {
                    m_buffer[row][column] = mBoard[row][column];
                }
            }
            for (UnsignedByte row = 0; row < m_dimension; row += 1) {
                delete[] mBoard[row];
            }
            delete[] mBoard;
            return maskId;
        }

        UnsignedByte numMasks = isMicro ? 4 : 8;
        UnsignedByte** maskedBoard[numMasks];
        size_t minScore = 0;
        UnsignedByte minId = 0;
        size_t maxScore = 0;
        UnsignedByte maxId = 0;
        for (UnsignedByte index = 0; index < numMasks; index += 1) {
            UnsignedByte mId = isMicro ? microMaskIdMap[index] : index;
            maskedBoard[index] = Mask(mId);
            size_t score = isMicro ?
                EvaluateMicro(maskedBoard[index]) :
                EvaluateCondition1(maskedBoard[index]) +
                EvaluateCondition2(maskedBoard[index]) +
                EvaluateCondition3(maskedBoard[index]) +
                EvaluateCondition4(maskedBoard[index]);
            if (minScore == 0 || minScore > score) {
                minScore = score;
                minId = index;
            }
            if (maxScore < score) {
                maxScore = score;
                maxId = index;
            }
        }

        UnsignedByte lasId = isMicro ? maxId : minId;
        for (UnsignedByte row = 0; row < m_dimension; row += 1) {
            for (UnsignedByte column = 0; column < m_dimension; column += 1) {
                m_buffer[row][column] = maskedBoard[lasId][row][column];
            }
        }

        for (UnsignedByte index = 0; index < numMasks; index += 1) {
            for (UnsignedByte row = 0; row < m_dimension; row += 1) {
                delete[] maskedBoard[index][row];
            }
            delete[] maskedBoard[index];
        }

        return lasId;
    }

    void GetFormatBits(ErrorCorrectionLevel level, UnsignedByte maskId, UnsignedByte* buffer)
    {
        static Unsigned2Bytes typeFormats[] = {
            0b1010100000100100, 0b1010001001001010, 0b1011110011111000, 0b1011011010010110, 0b1000101111110010, 0b1000000110011100, 0b1001111100101110, 0b1001010101000000,
            0b1110111110001000, 0b1110010111100110, 0b1111101101010100, 0b1111000100111010, 0b1100110001011110, 0b1100011000110000, 0b1101100010000010, 0b1101001011101100,
            0b0010110100010010, 0b0010011101111100, 0b0011100111001110, 0b0011001110100000, 0b0000111011000100, 0b0000010010101010, 0b0001101000011000, 0b0001000001110110,
            0b0110101010111110, 0b0110000011010000, 0b0111111001100010, 0b0111010000001100, 0b0100100101101000, 0b0100001100000110, 0b0101110110110100, 0b0101011111011010
        };

        UnsignedByte index = static_cast<UnsignedByte>(level) << 3;
        index |= maskId;

        Unsigned2Bytes value = typeFormats[index];
        UnsignedByte* valuePtr = (UnsignedByte*)&value;
        if constexpr (common::isLittleEndian) {
            buffer[0] = valuePtr[1];
            buffer[1] = valuePtr[0];
        } else {
            buffer[0] = valuePtr[0];
            buffer[1] = valuePtr[1];
        }
    }

    void GetMicroFormatBits(ErrorCorrectionLevel level, UnsignedByte version, UnsignedByte maskId, UnsignedByte* buffer)
    {
        static Unsigned2Bytes typeFormats[] {
            0b1000100010001010, 0b1000001011100100, 0b1001110001010110, 0b1001011000111000, 0b1010101101011100, 0b1010000100110010, 0b1011111110000000, 0b1011010111101110,
            0b1100111100100110, 0b1100010101001000, 0b1101101111111010, 0b1101000110010100, 0b1110110011110000, 0b1110011010011110, 0b1111100000101100, 0b1111001001000010,
            0b0000110110111100, 0b0000011111010010, 0b0001100101100000, 0b0001001100001110, 0b0010111001101010, 0b0010010000000100, 0b0011101010110110, 0b0011000011011000,
            0b0100101000010000, 0b0100000001111110, 0b0101111011001100, 0b0101010010100010, 0b0110100111000110, 0b0110001110101000, 0b0111110100011010, 0b0111011101110100
        };
        UnsignedByte index = GetMicroQRErrorCorrectionLevelValue(level, version);
        index = index << 2;
        index |= maskId;
        Unsigned2Bytes value = typeFormats[index];
        UnsignedByte* valuePtr = (UnsignedByte*)&value;
        if constexpr (common::isLittleEndian) {
            buffer[0] = valuePtr[1];
            buffer[1] = valuePtr[0];
        } else {
            buffer[0] = valuePtr[0];
            buffer[1] = valuePtr[1];
        }
    }

    void GetVersionBits(UnsignedByte version, UnsignedByte* buffer)
    {
        static Unsigned4Bytes versions[34] = {
            0b00011111001001010000000000000000,
            0b00100001011011110000000000000000,
            0b00100110101001100100000000000000,
            0b00101001001101001100000000000000,
            0b00101110111111011000000000000000,
            0b00110001110110001000000000000000,
            0b00110110000100011100000000000000,
            0b00111001100000110100000000000000,
            0b00111110010010100000000000000000,
            0b01000010110111100000000000000000,
            0b01000101000101110100000000000000,
            0b01001010100001011100000000000000,
            0b01001101010011001000000000000000,
            0b01010010011010011000000000000000,
            0b01010101101000001100000000000000,
            0b01011010001100100100000000000000,
            0b01011101111110110000000000000000,
            0b01100011101100010000000000000000,
            0b01100100011110000100000000000000,
            0b01101011111010101100000000000000,
            0b01101100001000111000000000000000,
            0b01110011000001101000000000000000,
            0b01110100110011111100000000000000,
            0b01111011010111010100000000000000,
            0b01111100100101000000000000000000,
            0b10000010011101010100000000000000,
            0b10000101101111000000000000000000,
            0b10001010001011101000000000000000,
            0b10001101111001111100000000000000,
            0b10010010110000101100000000000000,
            0b10010101000010111000000000000000,
            0b10011010100110010000000000000000,
            0b10011101010100000100000000000000,
            0b10100011000110100100000000000000,
        };

        Unsigned4Bytes value = versions[version - 7];
        UnsignedByte* valuePtr = (UnsignedByte*)&value;
        if constexpr (common::isLittleEndian) {
            buffer[0] = valuePtr[3];
            buffer[1] = valuePtr[2];
            buffer[2] = valuePtr[1];
        } else {
            buffer[0] = valuePtr[0];
            buffer[1] = valuePtr[1];
            buffer[2] = valuePtr[2];
        }
    }

    void SetReservedCell(UnsignedByte row, UnsignedByte column, UnsignedByte** buffer, UnsignedByte value)
    {
        UnsignedByte origin = buffer[row][column];
        UnsignedByte high = origin & 0xF0;
        if (high != BoardCell::FORMAT && high != BoardCell::VERSION) {
            ThrowException("Invalid reserved cell");
        }
        buffer[row][column] = high | (value & 0x0F);
    }

    void PlaceMicroFormat(UnsignedByte maskId, ErrorCorrectionInfo ecInfo)
    {
        UnsignedByte formatBits[2]; // Use first 15bits only
        GetMicroFormatBits(ecInfo.level, ecInfo.version, maskId, formatBits);

        for (UnsignedByte index = 0; index < 15; index += 1) {
            UnsignedByte byteIndex = index < 8 ? 0 : 1;
            UnsignedByte bitIndex = byteIndex == 0 ? index : index - 8;
            UnsignedByte mask = 0b10000000 >> bitIndex;
            UnsignedByte curByte = formatBits[byteIndex] & mask;
            UnsignedByte cell = curByte > 0 ? BoardCell::SET : BoardCell::UNSET;
            if (index < 8) {
                SetReservedCell(8, index + 1, m_buffer, cell);
            } else {
                SetReservedCell(15 - index, 8, m_buffer, cell);
            }
        }
    }

    void PlaceFormatAndVersion(UnsignedByte maskId, ErrorCorrectionInfo ecInfo)
    {
        UnsignedByte formatBits[2]; // Use first 15bits only
        GetFormatBits(ecInfo.level, maskId, formatBits);
        for (UnsignedByte index = 0; index < 15; index += 1) {
            UnsignedByte byteIndex = index < 8 ? 0 : 1;
            UnsignedByte bitIndex = byteIndex == 0 ? index : index - 8;
            UnsignedByte mask = 0b10000000 >> bitIndex;
            UnsignedByte curByte = formatBits[byteIndex] & mask;
            UnsignedByte cell = curByte > 0 ? BoardCell::SET : BoardCell::UNSET;
            if (index < 6) {
                SetReservedCell(8, index, m_buffer, cell);
                SetReservedCell(m_dimension - index - 1, 8, m_buffer, cell);
            } else if (index > 8) {
                SetReservedCell(8, m_dimension - (15 - index), m_buffer, cell);
                SetReservedCell(14 - index, 8, m_buffer, cell);
            } else {
                switch (index) {
                case 6:
                    SetReservedCell(8, index + 1, m_buffer, cell);
                    SetReservedCell(m_dimension - index - 1, 8, m_buffer, cell);
                    break;
                case 7:
                    SetReservedCell(8, index + 1, m_buffer, cell);
                    SetReservedCell(8, m_dimension - (15 - index), m_buffer, cell);
                    break;
                case 8:
                    SetReservedCell(14 - index + 1, 8, m_buffer, cell);
                    SetReservedCell(8, m_dimension - (15 - index), m_buffer, cell);
                    break;
                default:
                    break;
                }
            }
        }

        if (ecInfo.version < 7) return;

        UnsignedByte version[3];
        GetVersionBits(ecInfo.version, version);
        UnsignedByte row1 = m_dimension - 9;
        UnsignedByte col1 = 5;
        UnsignedByte row2 = 5;
        UnsignedByte col2 = m_dimension - 9;

        for (UnsignedByte index = 0; index < 18; index += 1) {
            UnsignedByte byteIndex = index / 8;
            UnsignedByte bitIndex = index % 8;
            UnsignedByte mask = 0b10000000 >> bitIndex;
            UnsignedByte curByte = version[byteIndex] & mask;
            UnsignedByte cell = curByte > 0 ? BoardCell::SET : BoardCell::UNSET;

            SetReservedCell(row1, col1, m_buffer, cell);
            if (row1 == (m_dimension - 11)) {
                row1 = m_dimension - 9;
                col1 -= 1;
            } else {
                row1 -= 1;
            }

            SetReservedCell(row2, col2, m_buffer, cell);
            if (col2 == (m_dimension - 11)) {
                col2 = m_dimension - 9;
                row2 -= 1;
            } else {
                col2 -= 1;
            }
        }
    }

private:
    UnsignedByte m_dimension{0};
    UnsignedByte ** m_buffer{nullptr};
};

struct QRMatrixStructuredAppend //todo, refactor
{
    /// Data segments
    QRMatrixSegment * segments{nullptr};
    /// Count of segments
    size_t count;
    /// Error correction info
    ErrorCorrectionLevel level;
    /// Optional. Limit minimum version
    /// (result version = max(minimum version, required version to fit data).
    UnsignedByte minVersion;
    /// Optional. Force to use given mask (0-7).
    /// Almost for test, you can ignore this.
    UnsignedByte maskId;
    /// Extra mode (MicroQR will be ignored. Not sure about FNC1.)
    QRMatrixExtraMode extraMode;

    QRMatrixStructuredAppend(QRMatrixSegment* segs, size_t segCount, ErrorCorrectionLevel ecLevel)
    {
        count = segCount;
        segments = new QRMatrixSegment[count];
        for (size_t index = 0; index < count; index++)
            segments[index] = segs[index];
        level = ecLevel;
        minVersion = 0;
        maskId = 0xFF;
        extraMode = QRMatrixExtraMode();
    }

    QRMatrixStructuredAppend()
    {
        count = 0;
        segments = new QRMatrixSegment[count];
        level = ErrorCorrectionLevel::LOW;
        minVersion = 0;
        maskId = 0xFF;
        extraMode = QRMatrixExtraMode();
    }

    QRMatrixStructuredAppend(QRMatrixStructuredAppend & other)
    {
        count = other.count;
        segments = new QRMatrixSegment [count];
        for (size_t index = 0; index < count; index++)
            segments[index] = other.segments[index];
        level = other.level;
        minVersion = other.minVersion;
        maskId = other.maskId;
        extraMode = other.extraMode;
    }

    ~QRMatrixStructuredAppend()
    {
        if (nullptr != segments) {
            delete[] segments;
            segments = nullptr;
        }
    }
    void operator=(QRMatrixStructuredAppend other)
    {
        count = other.count;
        segments = new QRMatrixSegment [count];
        for (size_t index = 0; index < count; index++)
            segments[index] = other.segments[index];
        level = other.level;
        minVersion = other.minVersion;
        maskId = other.maskId;
        extraMode = other.extraMode;
    }
};

namespace encoder {

inline static constexpr size_t ALPHA_NUM_MULTIPLICATION = 45;
inline static constexpr size_t ALPHA_NUM_PAIR_CHARS_BITS_LEN = 11;
inline static constexpr size_t ALPHA_NUM_SINGLE_CHAR_BITS_LEN = 6;

class AlphaNumericEncoder
{
public:
    /// Index of character in QR AlphaNumeric table.
    /// @ref: C++ std::string::find.
    static size_t GetIndexOfCharacter(unsigned char character)
    {
        static const char *table = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ $%*+-./:";
        const char* targetPtr = std::strchr(table, character);
        if (targetPtr != NULL) {
            return (size_t)(targetPtr - table);
        }
        return -1;
    }

    /// Calculate encoded value for pair of characters.
    /// @param pair length must be 1 or 2.
    static size_t EncodedValueOfPair(const UnsignedByte* pair, size_t length)
    {
        if (length < 1 || length > 2) {
            ThrowException("Pair must have 1 or 2 characters.");
        }
        size_t value1 = GetIndexOfCharacter(pair[0]);
        if (length == 1) return value1;

        size_t value2 = GetIndexOfCharacter(pair[1]);
        return value1 * ALPHA_NUM_MULTIPLICATION + value2;
    }

    /// Encode text.
    /// @return Number of written bits.
    static size_t Encode(const UnsignedByte* text, size_t length, UnsignedByte* buffer, size_t startIndex)
    {
        size_t index = 0;
        size_t bitIndex = startIndex;
        size_t charCount = 2;
        UnsignedByte pair[2];
        while (index < length) {
            pair[0] = text[index];
            size_t offset  = length - index;
            size_t encodedLen = 0;
            if (offset > 1) {
                charCount = 2;
                pair[1] = text[index + 1];
                encodedLen = ALPHA_NUM_PAIR_CHARS_BITS_LEN;
            } else {
                charCount = 1;
                encodedLen = ALPHA_NUM_SINGLE_CHAR_BITS_LEN;
            }
            index += charCount;
            size_t encodedData = EncodedValueOfPair(pair, charCount);
            size_t encodedSize = sizeof(size_t);
            size_t encodedOffset = encodedSize * 8 - encodedLen;
            CopyBits(
                (unsigned char*)&encodedData,
                encodedSize,
                encodedOffset,
                common::isLittleEndian,
                buffer,
                bitIndex,
                encodedLen
            );
            bitIndex += encodedLen;
        }
        return bitIndex - startIndex;
    }
};

class NumericEncoder
{
public:
    /// Encode text.
    /// @return Number of written bits.
    static size_t Encode(const UnsignedByte* text, size_t length, UnsignedByte* buffer, size_t startIndex)
    {
        int index = 0;
        size_t bitIndex = startIndex;
        while (index < length) {
            size_t groupLen;
            int offset  = length - index;
            if (offset > 2) {
                groupLen = 3;
            } else if (offset > 1) {
                groupLen = 2;
            } else {
                groupLen = 1;
            }
            std::string group = std::string((char*)text, index, groupLen);
            index += groupLen;
            size_t value = std::stoi(group);
            size_t bitLen;
            switch (groupLen) {
            case 3:
                bitLen = NUM_TRIPLE_DIGITS_BITS_LEN;
                break;
            case 2:
                bitLen = NUM_DOUBLE_DIGITS_BITS_LEN;
                break;
            case 1:
                bitLen = NUM_SINGLE_DIGIT_BITS_LEN;
                break;
            default:
                bitLen = 0;
                break;
            }
            size_t encodedSize = sizeof(size_t);
            size_t encodedOffset = encodedSize * 8 - bitLen;
            CopyBits(
                (UnsignedByte*)&value,
                encodedSize,
                encodedOffset,
                common::isLittleEndian,
                buffer,
                bitIndex,
                bitLen
            );
            bitIndex += bitLen;
        }
        return bitIndex - startIndex;
    }
};

class KanjiEncoder
{
public:
    /// Encode text.
    /// @return Number of written bits.
    static size_t Encode(const UnsignedByte* text, size_t length, UnsignedByte* buffer, size_t startIndex)
    {
        size_t bitIndex = startIndex;
        for (size_t index = 0; index < length; index += 2) {
            UnsignedByte* ptr = (UnsignedByte*)text + index;
            Unsigned2Bytes charWord = 0;
            UnsignedByte* charWordPtr = (UnsignedByte*)&charWord;
            if constexpr (common::isLittleEndian) {
                *charWordPtr = *(ptr + 1);
                *(charWordPtr + 1) = *ptr;
            } else {
                *charWordPtr = *ptr;
                *(charWordPtr + 1) = *(ptr + 1);
            }
            Unsigned2Bytes offset = 0;
            if (charWord >= 0x8140 && charWord <= 0x9FFC) {
                offset = 0x8140;
            } else if (charWord >= 0xE040 && charWord <= 0xEBBF) {
                offset = 0xC140;
            } else {
                ThrowException("Unknown or unsupported kanji character.");
            }
            charWord = charWord - offset;
            Unsigned2Bytes lsByte = 0;
            Unsigned2Bytes msByte = 0;
            if constexpr (common::isLittleEndian) {
                lsByte = (Unsigned2Bytes)*charWordPtr;
                msByte = (Unsigned2Bytes)*(charWordPtr + 1);
            } else {
                lsByte = (Unsigned2Bytes)*(charWordPtr + 1);
                msByte = (Unsigned2Bytes)*charWordPtr;
            }
            charWord = (msByte * 0xC0) + lsByte;
            CopyBits(charWordPtr, 2, 3 /* 16 - 13 */, common::isLittleEndian, buffer, bitIndex, 13);
            bitIndex += 13;
        }
        return bitIndex - startIndex;
    }
};

} // namespace encoder

class QRMatrixEncoder
{
public:

    /// Get QR Version (dimension) to encode given data.
    /// @return 0 if no suiversion
    static UnsignedByte GetVersion(
        /// Array of segments to be encoded
        QRMatrixSegment* segments,
        /// Number of segments
        size_t count,
        /// Error correction info
        ErrorCorrectionLevel level,
        /// Extra mode
        QRMatrixExtraMode extraMode = QRMatrixExtraMode(),
        /// Is this symbol a part of Structured Append
        bool isStructuredAppend = false
    )
    {
        size_t segCount = 0;
        for (size_t index = 0; index < count; index += 1) {
            if (segments[index].Length() > 0) {
                segCount += 1;
            }
        }
        if (segCount == 0) {
            return 0;
        }
        try {
            ErrorCorrectionInfo ecInfo = FindVersion(segments, count, level, 0, extraMode, isStructuredAppend);
            if (ecInfo.version == 0) {
                return 0;
            }
            return ecInfo.version;
        }
        catch (const Exception & exception){
            return 0;
        }
    }

    /// Encode single QR symbol
    static QRMatrixBoard Encode(
        /// Array of segments to be encoded
        QRMatrixSegment* segments,
        /// Number of segments
        size_t count,
        /// Error correction info
        ErrorCorrectionLevel level,
        /// Extra mode
        QRMatrixExtraMode extraMode = QRMatrixExtraMode(),
        /// Optional. Limit minimum version
        /// (result version = max(minimum version, required version to fit data).
        UnsignedByte minVersion = 0,
        /// Optional. Force to use given mask (0-7).
        /// Almost for test, you can ignore this.
        UnsignedByte maskId = 0xFF
    )
    {
        return EncodeSingle(segments, count, level, extraMode, minVersion, maskId, UnsignedByte{0}, UnsignedByte{0}, UnsignedByte{0});
    }

    /// Encode Structured Append QR symbols
    /// @return Array of QRMatrixBoard (should be deleted when done).
    static QRMatrixBoard* Encode(
        /// Array of data parts to be encoded
        QRMatrixStructuredAppend* parts,
        /// Number of parts
        size_t count
    )
    {
        if (count > 16) {
            ThrowException("Structured Append only accepts 16 parts maximum");
        }
        if (count == 0) {
            ThrowException("No input.");
        }
    //    if (count == 1) {
            // Should change to encode single QR symbol or throw error?
            // But no rule prevents to make a Structured Append QR symbol single part
    //    }
        UnsignedByte parity = 0;
        for (UnsignedByte index = 0; index < count; index += 1) {
            QRMatrixStructuredAppend part = parts[index];
            for (size_t segIndex = 0; segIndex < part.count; segIndex += 1) {
                QRMatrixSegment segment = part.segments[segIndex];
                for (size_t idx = 0; idx < segment.Length(); idx += 1) {
                    if (index == 0 && segIndex == 0 && idx == 0) {
                        parity = segment.Data()[idx];
                    } else {
                        parity = parity ^ segment.Data()[idx];
                    }
                }
            }
        }
        QRMatrixBoard* result = new QRMatrixBoard[count];
        for (UnsignedByte index = 0; index < count; index += 1) {
            QRMatrixStructuredAppend part = parts[index];
            try {
                result[index] = EncodeSingle(
                    part.segments, part.count,
                    part.level, part.extraMode,
                    part.minVersion, part.maskId,
                    index, count, parity
                );
            } catch (const Exception & exception) {
                delete[] result;
                throw exception;
                return nullptr;
            }
        }
        return result;
    }


    // Calculate encoded data bits count,
    /// include Mode Indicator and ECI header bits,
    /// exclude Characters Count bits
    static size_t CalculateEncodedDataBitsCount(QRMatrixSegment* segments, size_t count)
    {
        size_t totalDataBitsCount = 0;
        for (size_t index = 0; index < count; index += 1) {
            QRMatrixSegment segment = segments[index];
            if (segment.Length() == 0) {
                continue;
            }
            // Mode indicator
            totalDataBitsCount += 4;
            // ECI
            if (segment.isEciHeaderRequired()) {
                totalDataBitsCount += 4; // ECI Header
                if (segment.ECI() <= 128) {
                    totalDataBitsCount += 8; // 1 byte ECI Indicator
                } else if (segment.ECI() <= 16383) {
                    totalDataBitsCount += 16; // 2 bytes ECI Indicator
                } else {
                    totalDataBitsCount += 24; // 3 bytes ECI Indicator
                }
            }
            // Data
            switch (segment.Mode()) {
            case EncodingMode::NUMERIC: {
                // 3 characters encoded in 10 bits (each character is 1 byte)
                size_t numberOfGroups = (segment.Length() / 3);
                totalDataBitsCount += numberOfGroups * NUM_TRIPLE_DIGITS_BITS_LEN;
                // Remaining chars
                UnsignedByte remainChars = segment.Length() % 3;
                switch (remainChars) {
                case 1:
                    totalDataBitsCount += NUM_SINGLE_DIGIT_BITS_LEN;
                    break;
                case 2:
                    totalDataBitsCount += NUM_DOUBLE_DIGITS_BITS_LEN;
                default:
                    break;
                }
            }
            break;
            case EncodingMode::ALPHA_NUMERIC: {
                // 2 characters encoded in 11 bits (each character is 1 byte)
                // Remaining character encoded in 6 bits.
                size_t numberOfGroups = segment.Length() / 2;
                size_t remaining = segment.Length() % 2;
                totalDataBitsCount += ALPHA_NUM_PAIR_CHARS_BITS_LEN * numberOfGroups + ALPHA_NUM_SINGLE_CHAR_BITS_LEN * remaining;
            }
            break;
            case EncodingMode::KANJI:
                // 2 bytes per Kanji character.
                // Each character is encoded in 13 bits.
                totalDataBitsCount += (segment.Length() / 2) * 13;
                break;
            case EncodingMode::BYTE:
                totalDataBitsCount += segment.Length() * 8;
                break;
            }
        }
        return totalDataBitsCount;
    }

    /// Find QR Version & its properties
    static ErrorCorrectionInfo  FindVersion(QRMatrixSegment* segments, size_t count, ErrorCorrectionLevel level, UnsignedByte minVersion, QRMatrixExtraMode extraMode, bool isStructuredAppend)
    {
        // This is total estimated bits of data, excluding bits for Characters Count,
        // because each version requires difference number of Characters Count bits.
        size_t totalDataBitsCount = CalculateEncodedDataBitsCount(segments, count);
        if (isStructuredAppend) {
            totalDataBitsCount += 20; // 4 bits header, 4 bits position, 4 bits total number, 1 byte parity
        }
        if (extraMode.mode == EncodingExtraMode::FNC1_FIRST) {
            totalDataBitsCount += 4; // 4 bits FNC1 indicator
        } else  if (extraMode.mode == EncodingExtraMode::FNC1_SECOND) {
            totalDataBitsCount += 12; // 4 bits FNC1 indicator, 8 bits Application Indicator
        }
        // Go for each version, get total bits and check
        bool isMicro = (extraMode.mode == EncodingExtraMode::MICRO_QR && !isStructuredAppend);
        UnsignedByte maxVer = isMicro ? MICROQR_MAX_VERSION : QR_MAX_VERSION;
        for (UnsignedByte version = (minVersion > 0 && minVersion <= maxVer) ? minVersion : 1; version <= maxVer; version += 1) {
            if (version > 1 && isMicro && level == ErrorCorrectionLevel::HIGH) {
                ThrowException("Error Correction Level High not available in MicroQR.");
            }
            ErrorCorrectionInfo info = isMicro ?
                ErrorCorrectionInfo::GetMicroErrorCorrectionInfo(version, level) :
                ErrorCorrectionInfo::GetErrorCorrectionInfo(version, level);
            size_t capacity = info.codewords * 8;
            if (capacity == 0) {
                break;
            }
            if (totalDataBitsCount >= capacity) {
                // Not need to count total bits, just go to next version
                continue;
            }
            // Calculate total number of bits for Characters Count for this version
            size_t totalBits = totalDataBitsCount;
            bool hasAlpha = false;
            bool hasKanji = false;
            bool hasByte = false;
            for (size_t index = 0; index < count; index += 1) {
                QRMatrixSegment segment = segments[index];
                if (segment.Length() == 0) {
                    continue;
                }
                switch (segment.Mode()) {
                case EncodingMode::ALPHA_NUMERIC:
                    hasAlpha = true;
                    break;
                case EncodingMode::BYTE:
                    hasByte = true;
                    break;
                case EncodingMode::KANJI:
                    hasKanji = true;
                    break;
                default:
                    break;
                }
                totalBits += isMicro ?
                    GetMicroCharactersCountIndicatorLength(version, segment.Mode()) :
                    GetCharactersCountIndicatorLength(version, segment.Mode());
            }
            // Check if this QR version bits capacity is enough for required total bits
            if (totalBits <= capacity) {
                UnsignedByte finalVer = version;
                if (isMicro) {
                    if (finalVer < 2 && hasAlpha) {
                        // AlphaNumeric Mode is not available with M1
                        finalVer = 2;
                    }
                    if (finalVer < 3 && (hasByte || hasKanji)) {
                        // Byte/Kanji Mode is not available with <= M2
                        finalVer = 3;
                    }
                }
                if (finalVer == version) {
                    return info;
                }
                return isMicro ?
                    ErrorCorrectionInfo::GetMicroErrorCorrectionInfo(finalVer, level) :
                    ErrorCorrectionInfo::GetErrorCorrectionInfo(finalVer, level);
            }
        }
        return ErrorCorrectionInfo();
    }

    /// Encode ECI Indicator
    static UnsignedByte* EncodeEciIndicator(Unsigned4Bytes indicator, UnsignedByte* resultLength)
    {
        UnsignedByte* indicatorPtr = (UnsignedByte*)&indicator;
        UnsignedByte* result;
        if (indicator <= 127) {
            *resultLength = 1;
            result = new UnsignedByte [1];
            if constexpr (common::isLittleEndian) {
                result[0] = indicatorPtr[0] & 0x7F;
            } else {
                result[0] = indicatorPtr[3] & 0x7F;
            }
        } else if (indicator <= 16383) {
            *resultLength = 2;
            result = new UnsignedByte [2];
            if constexpr (common::isLittleEndian) {
                result[0] = (indicatorPtr[1] & 0x3F) | 0x80;
                result[1] = indicatorPtr[0];
            } else {
                result[0] = (indicatorPtr[2] & 0x3F) | 0x80;
                result[1] = indicatorPtr[3];
            }
        } else if (indicator <= 999999) {
            *resultLength = 3;
            result = new UnsignedByte [3];
            if constexpr (common::isLittleEndian) {
                result[0] = (indicatorPtr[2] & 0x1F) | 0xC0;
                result[1] = indicatorPtr[1];
                result[2] = indicatorPtr[0];
            } else {
                result[0] = (indicatorPtr[1] & 0x1F) | 0xC0;
                result[1] = indicatorPtr[2];
                result[2] = indicatorPtr[3];
            }
        } else {
            ThrowException("Invalid ECI Indicator");
        }
        return result;
    }

    static QRMatrixBoard EncodeSingle(
        QRMatrixSegment* segments,
        size_t count,
        ErrorCorrectionLevel level,
        QRMatrixExtraMode extraMode,
        UnsignedByte minVersion,
        UnsignedByte maskId,
        UnsignedByte sequenceIndex,
        UnsignedByte sequenceTotal,
        UnsignedByte parity
    ) 
    {
        size_t segCount = 0;
        for (size_t index = 0; index < count; index += 1) {
            if (segments[index].Length() > 0) {
                segCount += 1;
            }
        }
        if (segCount == 0) {
            ThrowException("No input.");
        }
        if (extraMode.mode == EncodingExtraMode::FNC1_SECOND) {
            bool isValid = false;
            switch (extraMode.appIndicatorLength) {
            case 1: {
                UnsignedByte value = extraMode.appIndicator[0];
                isValid = (value >= 'a' && value <= 'z') || (value >= 'A' && value <= 'Z');
            }
                break;
            case 2:
            {
                UnsignedByte value1 = extraMode.appIndicator[0];
                UnsignedByte value2 = extraMode.appIndicator[1];
                isValid = (value1 >= '0' && value1 <= '9') && (value2 >= '0' && value2 <= '9');
            }
                break;
            default:
                break;
            }
            if (not isValid) {
                ThrowException("Invalid Application Indicator for FNC1 Second Position mode");
            }
        }
        bool isStructuredAppend = sequenceTotal > 0 && sequenceTotal <= 16;
        ErrorCorrectionInfo ecInfo = FindVersion(segments, count, level, minVersion, extraMode, isStructuredAppend);
        if (ecInfo.version == 0) {
            ThrowException("Unable to find suitable QR version.");
        }
        if (isStructuredAppend && extraMode.mode == EncodingExtraMode::MICRO_QR) {
            extraMode = QRMatrixExtraMode();
        }
        // Allocate
        UnsignedByte* buffer = Allocate(ecInfo.codewords);
        size_t bitIndex = 0;
        // Structured append
        if (isStructuredAppend) {
            UnsignedByte bufByte = 0b0011;
            CopyBits(&bufByte, 1, 4, false, buffer, bitIndex, 4);
            bitIndex += 4;
            CopyBits(&sequenceIndex, 1, 4, false, buffer, bitIndex, 4);
            bitIndex += 4;
            bufByte = sequenceTotal - 1;
            CopyBits(&bufByte, 1, 4, false, buffer, bitIndex, 4);
            bitIndex += 4;
            CopyBits(&parity, 1, 0, false, buffer, bitIndex, 8);
            bitIndex += 8;
        }
        // Encode data
        for (size_t index = 0; index < count; index += 1) {
            EncodeSegment(buffer, segments[index], index, level, ecInfo, &bitIndex, extraMode);
        }
        // Finish
        return FinishEncodingData(buffer, ecInfo, &bitIndex, maskId, extraMode);
    }

    static QRMatrixBoard FinishEncodingData(
        UnsignedByte* buffer,
        ErrorCorrectionInfo ecInfo,
        size_t* bitIndex,
        UnsignedByte maskId,
        QRMatrixExtraMode extraMode
    )
    {
        bool isMicro = (extraMode.mode == EncodingExtraMode::MICRO_QR);
        bool isMicroV13 = isMicro && ((ecInfo.version == 1) || ecInfo.version == 3);
        size_t bufferLen = ecInfo.codewords;
        size_t bufferBitsLen = bufferLen * 8;
        if (isMicroV13) {
            bufferBitsLen -= 4;
        }
        /// Terminator
        size_t terminatorLength = isMicro ? GetMicroTerminatorLength(ecInfo.version) : 4;
        for (size_t index = 0; *bitIndex < bufferBitsLen && index < terminatorLength; index += 1) {
            *bitIndex += 1;
        }

        /// Make data multiple by 8
        if (isMicroV13) {
            while ((*bitIndex < bufferBitsLen - 4) && (*bitIndex % 8 != 0)) {
                *bitIndex += 1;
            }
        } else {
            while (*bitIndex < bufferBitsLen && *bitIndex % 8 != 0) {
                *bitIndex += 1;
            }
        }

        /// Fill up
        UnsignedByte byteFilling1 = 0b11101100; // 0xEC
        UnsignedByte byteFilling2 = 0b00010001; // 0x11
        UnsignedByte curByteFilling = byteFilling1;
        while (*bitIndex < bufferBitsLen - (isMicroV13 ? 4 : 0)) {
            CopyBits(
                (UnsignedByte*)&curByteFilling,
                1, // size of curByteFilling
                0, // Bit index of curByteFilling
                false, // copy 1 byte source, so we don't need to care byte order here
                buffer, // destination
                *bitIndex, // destination bit index
                8
            );
            *bitIndex += 8;
            if (curByteFilling == byteFilling1) {
                curByteFilling = byteFilling2;
            } else {
                curByteFilling = byteFilling1;
            }
        }
        *bitIndex = bufferBitsLen;

        // Error corrections
        UnsignedByte** ecBuffer = GenerateErrorCorrections(buffer, ecInfo);

        // Interleave or not
        if (ecInfo.EcBlockTotalCount() > 1) {
            UnsignedByte *interleave = Interleave(buffer, ecInfo);
            UnsignedByte *ecInterleave = Interleave(ecBuffer, ecInfo);

            QRMatrixBoard board(interleave, ecInterleave, ecInfo, maskId, isMicro);

            delete[] interleave;
            delete[] ecInterleave;
            Clean(buffer, ecBuffer, ecInfo);
            return board;
        } else {
            QRMatrixBoard board(buffer, ecBuffer[0], ecInfo, maskId, isMicro);
            Clean(buffer, ecBuffer, ecInfo);
            return board;
        }
    }

    /// Generate Error correction bytes
    /// @return Binary data (should be deleted on unused)
    static UnsignedByte* GenerateErrorCorrections(
        /// Bytes from previous (encode data) step
        UnsignedByte* encodedData,
        /// EC Info from previous step
        ErrorCorrectionInfo ecInfo,
        /// Group number: 0, 1
        UnsignedByte group,
        /// Block number: 0, ...
        UnsignedByte block
    )
    {
        UnsignedByte maxGroup = ecInfo.GroupCount();
        if (group >= maxGroup) {
            ThrowException("Invalid group number");
        }
        UnsignedByte maxBlock = 0;
        size_t offset = 0;
        size_t blockSize = 0;
        switch (group) {
        case 0:
            maxBlock = ecInfo.group1Blocks;
            blockSize = ecInfo.group1BlockCodewords;
            offset = block * ecInfo.group1BlockCodewords;
            break;
        case 1:
            maxBlock = ecInfo.group2Blocks;
            blockSize = ecInfo.group2BlockCodewords;
            offset = ecInfo.group1Blocks * ecInfo.group1BlockCodewords + ecInfo.group2BlockCodewords * block;
            break;
        default:
            break;
        }
        if (block >= maxBlock) {
            ThrowException("Invalid block number");
        }

        Polynomial mesg(blockSize);
        for (size_t index = 0; index < blockSize; index += 1) {
            mesg.terms[index] = encodedData[offset + index];
        }
        Polynomial ecc = mesg.GetErrorCorrections(ecInfo.ecCodewordsPerBlock);
        UnsignedByte* result = Allocate(ecc.length);
        for (size_t index = 0; index < ecc.length; index += 1) {
            result[index] = ecc.terms[index];
        }
        return result;
    }

    /// Generate Error correction bytes
    /// @return Binary data (should be deleted on unused). 2D array contains EC bytes for all blocks.
    static UnsignedByte** GenerateErrorCorrections(
        /// Bytes from previous (encode data) step
        UnsignedByte* encodedData,
        /// EC Info from previous step
        ErrorCorrectionInfo ecInfo
    ) {
        UnsignedByte maxGroup = ecInfo.GroupCount();
        UnsignedByte** result = new UnsignedByte* [ecInfo.EcBlockTotalCount()];
        size_t blockIndex = 0;
        for (UnsignedByte group = 0; group < maxGroup; group += 1) {
            UnsignedByte maxBlock = 0;
            switch (group) {
            case 0:
                maxBlock = ecInfo.group1Blocks;
                break;
            case 1:
                maxBlock = ecInfo.group2Blocks;
                break;
            default:
                break;
            }
            for (UnsignedByte block = 0; block < maxBlock; block += 1) {
                result[blockIndex] = GenerateErrorCorrections(encodedData, ecInfo, group, block);
                blockIndex += 1;
            }
        }
        return result;
    }

    /// Interleave data codeworks
    /// Throw error if QR has only 1 block in total (check ecInfo before call this function)
    /// @return Binary data (should be deleted on unused)
    static UnsignedByte* Interleave(
        /// Encoded data
        UnsignedByte* encodedData,
        /// EC info
        ErrorCorrectionInfo ecInfo
    )
    {
        if (ecInfo.EcBlockTotalCount() == 1) {
            ThrowException("Interleave not required");
        }
        UnsignedByte* result = new UnsignedByte [ecInfo.codewords];
        size_t blockCount = ecInfo.EcBlockTotalCount();
        UnsignedByte* blockPtr[blockCount];
        UnsignedByte* blockEndPtr[blockCount];

        UnsignedByte groupCount = ecInfo.GroupCount();
        size_t blockIndex = 0;

        for (UnsignedByte group = 0; group < groupCount; group += 1) {
            UnsignedByte blockCount = ecInfo.group1Blocks;
            size_t blockSize = ecInfo.group1BlockCodewords;
            if (group == 1) {
                blockCount = ecInfo.group2Blocks;
                blockSize = ecInfo.group2BlockCodewords;
            }
            for (UnsignedByte block = 0; block < blockCount; block += 1) {
                size_t offset = 0;
                if (group == 1) {
                    offset += ecInfo.group1Blocks * ecInfo.group1BlockCodewords;
                }
                offset += block * blockSize;
                blockPtr[blockIndex] = &encodedData[offset];
                blockEndPtr[blockIndex] = &encodedData[offset + blockSize - 1];
                blockIndex += 1;
            }
        }

        blockIndex = 0;
        size_t resIndex = 0;
        bool found = true;
        while (found) {
            result[resIndex] = *blockPtr[blockIndex];
            resIndex += 1;
            blockPtr[blockIndex] += 1;

            size_t loopCount = 0;
            found = false;
            while (loopCount < blockCount) {
                blockIndex += 1;
                if (blockIndex >= blockCount) {
                    blockIndex = 0;
                }
                if (blockPtr[blockIndex] <= blockEndPtr[blockIndex]) {
                    found = true;
                    break;
                }
                loopCount += 1;
            }
        }
        return result;
    }

    /// Interleave error correction codeworks
    /// Throw error if QR has only 1 block in total (check ecInfo before call this function)
    /// @return Binary data (should be deleted on unused)
    static UnsignedByte* Interleave(
        /// Error correction data
        UnsignedByte** data,
        /// EC info
        ErrorCorrectionInfo ecInfo
    )
    {
        if (ecInfo.EcBlockTotalCount() == 1) {
            ThrowException("Interleave not required");
        }
        size_t blockCount = ecInfo.EcBlockTotalCount();
        UnsignedByte* result = new UnsignedByte[ecInfo.EcCodewordsTotalCount()];
        size_t resIndex = 0;
        for (size_t index = 0; index < ecInfo.ecCodewordsPerBlock; index += 1) {
            for (size_t jndex = 0; jndex < blockCount; jndex += 1) {
                result[resIndex] = data[jndex][index];
                resIndex += 1;
            }
        }
        return result;
    }

    static void Clean(UnsignedByte* buffer, UnsignedByte** ecBuffer, ErrorCorrectionInfo ecInfo)
    {
        for (size_t index = 0; index < (ecInfo.group1Blocks + ecInfo.group2Blocks); index += 1) {
            delete[] ecBuffer[index];
        }
        delete[] ecBuffer;
        delete[] buffer;
    }

    /// Encode segments into buffer
    static void EncodeSegment(
        UnsignedByte* buffer,
        QRMatrixSegment segment,
        size_t segmentIndex,
        ErrorCorrectionLevel level,
        ErrorCorrectionInfo ecInfo,
        size_t* bitIndex,
        QRMatrixExtraMode extraMode)
    {
        if (segment.Length() == 0) {
            return;
        }
        bool isMicro = extraMode.mode == EncodingExtraMode::MICRO_QR;
        size_t uintSize = sizeof(size_t);
        // ECI Header if enable
        if (!isMicro && segment.isEciHeaderRequired()) {
            UnsignedByte eciLen = 0;
            UnsignedByte * eciHeader = EncodeEciIndicator(segment.ECI(), &eciLen);
            // 4 bits of ECI mode indicator
            UnsignedByte eciModeHeader = 0b0111;
            CopyBits(&eciModeHeader, 1, 4, false, buffer, *bitIndex, 4);
            *bitIndex += 4;
            // ECI indicator
            CopyBits(
                eciHeader,
                eciLen, // size of eciUtf8Indicator
                0, // bit of eciUtf8Indicator to start copy
                false, // `eciHeader` is always normal order (big endian)
                buffer,
                *bitIndex,
                8 * eciLen // number of bits to be copied
                );
            *bitIndex += 8 * eciLen;
            delete[] eciHeader;
        }
        if (segmentIndex == 0) {
            if (extraMode.mode == EncodingExtraMode::FNC1_FIRST) {
                UnsignedByte fnc1Header = 0b0101;
                CopyBits(&fnc1Header, 1, 4, false, buffer, *bitIndex, 4);
                *bitIndex += 4;
            } else if (extraMode.mode == EncodingExtraMode::FNC1_SECOND) {
                UnsignedByte fnc1Header = 0b1001;
                CopyBits(&fnc1Header, 1, 4, false, buffer, *bitIndex, 4);
                *bitIndex += 4;
                fnc1Header = 0;
                switch (extraMode.appIndicatorLength) {
                case 1:
                    fnc1Header = extraMode.appIndicator[0] + 100;
                    break;
                case 2: {
                    std::string numStr = std::string((char*)extraMode.appIndicator, 2);
                    fnc1Header = std::stoi(numStr);
                }
                    break;
                default:
                    break;
                }
                CopyBits(&fnc1Header, 1, 0, false, buffer, *bitIndex, 8);
                *bitIndex += 8;
            }
        }
        // Bits of segment mode indicator
        UnsignedByte numberOfModeBits = isMicro ? GetMicroModeIndicatorLength(ecInfo.version, segment.Mode()) : 4;
        if (numberOfModeBits > 0) {
            UnsignedByte mode = isMicro ? GetMicroQREncodingModeValue(segment.Mode()) : static_cast<UnsignedByte>(segment.Mode());
            CopyBits(&mode, 1, 8 - numberOfModeBits, false, buffer, *bitIndex, numberOfModeBits);
            *bitIndex += numberOfModeBits;
        }
        // Character counts bits
        size_t charCountIndicatorLen = isMicro ?
            GetMicroCharactersCountIndicatorLength(ecInfo.version, segment.Mode()) :
            GetCharactersCountIndicatorLength(ecInfo.version, segment.Mode());
        size_t charCount = 0;
        switch (segment.Mode()) {
        case EncodingMode::NUMERIC:
        case EncodingMode::ALPHA_NUMERIC:
        case EncodingMode::BYTE:
            charCount = segment.Length();
            break;
        case EncodingMode::KANJI:
            charCount = segment.Length() / 2;
            break;
        }
        CopyBits(
            (UnsignedByte*)&charCount,
            uintSize,
            uintSize * 8 - charCountIndicatorLen,
            common::isLittleEndian,
            buffer,
            *bitIndex,
            charCountIndicatorLen
        );
        *bitIndex += charCountIndicatorLen;

        switch (segment.Mode()) {
        case EncodingMode::NUMERIC:
            *bitIndex += encoder::NumericEncoder::Encode(segment.Data(), segment.Length(), buffer, *bitIndex);
            break;
        case EncodingMode::ALPHA_NUMERIC:
            *bitIndex += encoder::AlphaNumericEncoder::Encode(segment.Data(), segment.Length(), buffer, *bitIndex);
            break;
        case EncodingMode::BYTE: {
            UnsignedByte* bytes = segment.Data();
            for (size_t idx = 0; idx < segment.Length(); idx += 1) {
                CopyBits(bytes, 1, 0, false, buffer, *bitIndex, 8);
                *bitIndex += 8;
                bytes += 1;
            }
        }
            break;
        case EncodingMode::KANJI:
            *bitIndex += encoder::KanjiEncoder::Encode(segment.Data(), segment.Length(), buffer, *bitIndex);
            break;
        }
    }
};

inline void FillCell(boost::gil::gray8_image_t::view_t & view, size_t dimension, size_t row, size_t column, UnsignedByte scale, UnsignedByte quiteZone)
{
    size_t rowOffset = (row + quiteZone) * scale * dimension;
    size_t colOffset = (column + quiteZone) * scale;
    for (size_t index = 0; index < scale; index += 1) {
        size_t pic = rowOffset + index * dimension + colOffset;
        for (unsigned jndex = 0; jndex < scale; jndex += 1) {
            // image[pic] = 0;
            pic += 1;
        }
    }

    //todo
}

inline void makeQR(const QRMatrixBoard & board, std::string_view path, bool isMicro)
{
    UnsignedByte quietZone = isMicro ? 2 : 4;
    static UnsignedByte scale = 10;
    size_t dimension = (board.Dimension() + 2 * quietZone) * scale;
    
    if (not fs::CreateDir(fs::DirName(path)))
        ThrowException("failed to create directory " + fs::DirName(path).string());
    
    boost::gil::gray8_image_t img(dimension, dimension);
    boost::gil::gray8_image_t::view_t v = boost::gil::view(img);
    boost::gil::fill_pixels(v, boost::gil::gray8_pixel_t(255));

    // Draw QR board (black cells)
    for (size_t row = 0; row < board.Dimension(); row++) {
        for (size_t column = 0; column < board.Dimension(); column++) {
            UnsignedByte cell = board.Buffer()[row][column];
            UnsignedByte low = cell & BoardCell::LOWMASK;
            if (low == BoardCell::SET)
                FillCell(v, dimension, row, column, scale, quietZone);
        }
    }
    
    boost::gil::write_view(path.data(), boost::gil::view(img), boost::gil::png_tag());
}

} // namespace detail

bool EncodeQR(std::string_view path, const UnsignedByte * raw, EncodingMode mode, size_t eci, bool isMicro, std::string * errMsg = nullptr)
{
    using namespace detail;
    try {
        ErrorCorrectionLevel level = isMicro ? ErrorCorrectionLevel::LOW : ErrorCorrectionLevel::HIGH;
        QRMatrixSegment segment(mode, raw, std::strlen((const char*)raw), eci != 0 ? eci : DEFAULT_ECI_ASSIGMENT_VALUE);
        QRMatrixSegment segments[] = {segment};
        QRMatrixExtraMode extra = isMicro ? QRMatrixExtraMode(EncodingExtraMode::MICRO_QR) : QRMatrixExtraMode();
        QRMatrixBoard board = QRMatrixEncoder::Encode(segments, 1, level, extra);
        std::cout << board.Description() << std::endl;
        makeQR(board, path, isMicro);
    } catch (const Exception & e) {
        if (errMsg) *errMsg = e.what();
        return false;
    }
    return true;
}

} // namespace generic::img::qr