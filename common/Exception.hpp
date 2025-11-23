/**
 * @file Exception.hpp
 * @author bwu
 * @brief Exception handle
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "Macros.hpp"
#include <system_error>
#include <exception>
#include <string>

#ifdef GENERIC_NO_EXCEPTION
#   define GENERIC_TRY
#   define GENERIC_CATCH
#   define GENERIC_THROW(ex) do{} while(0);
#else
#   define GENERIC_TRY try
#   define GENERIC_CATCH catch (const std::exception & ) {}
#   define GENERIC_THROW(ex) throw(ex);
#endif

#ifdef GENERIC_NO_ASSERTION
#   define GENERIC_ASSERT(ex) do{} while(0)
#   define GENERIC_ASSERT_MSG(ex, msg)                   \
    do {                                                 \
        static_cast<void>(sizeof(ex));                   \
        static_cast<void>(sizeof(msg));                  \
    } while (0)
#else
#   undef NDEBUG
#   include <cassert>
#   include <cstdio>
#   define GENERIC_ASSERT(ex) assert(ex)
#   define GENERIC_ASSERT_MSG(ex, msg)                    \
    if (not (ex)) {                                       \
        fprintf(stderr, "%s:%d ", __FILE__, __LINE__);    \
        fprintf(stderr, "Assertion '%s' failed: ", #ex);  \
        fprintf(stderr, "%s", msg);                       \
        fprintf(stderr, "!\n");                           \
        std::abort();                                     \
    }   
#endif

///@brief generic library namespace
namespace generic{

class Exception : public std::exception
{
public:
/**
 * @brief Brief description of Exception.
 * @param m_msg(std::move(msg)
 * @return explicit
 */
    explicit Exception(std::string msg) : m_msg(std::move(msg)) {}
    Exception(std::string_view msg, int errCode)
    {
        auto ec = std::error_code(errCode, std::generic_category());
        m_msg.append(std::system_error(ec, msg.data()).what());
    }
    const char * what() const noexcept override { return m_msg.c_str(); }

private:
    std::string m_msg;
};

inline void ThrowException(std::string msg)
{
    GENERIC_THROW(Exception(std::move(msg)));
}

inline void ThrowException(std::string_view msg, int errCode)
{
    GENERIC_THROW(Exception(msg, errCode));
}

}//namespace generic
