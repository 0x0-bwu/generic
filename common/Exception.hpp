/**
 * @file Exception.hpp
 * @author bwu
 * @brief Exception handle
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_COMMON_EXCEPTION_HPP
#define GENERIC_COMMON_EXCEPTION_HPP
#include "Macros.hpp"
#include <system_error>
#include <exception>
#include <cassert>
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

#ifdef GENERIC_NO_ASSERT
#   define GENERIC_ASSERT(ex) do{} while(0);
#else
#   define GENERIC_ASSERT(ex) assert(ex);
#endif

///@brief generic library namespace
namespace generic{
///@brief common classes and functions
namespace common {

class Exception : public std::exception
{
public:
    explicit Exception(std::string msg) : m_msg(std::move(msg)) {}
    Exception(const std::string & msg, int errCode)
    {
        auto ec = std::error_code(errCode, std::generic_category());
        m_msg.append(std::system_error(ec, msg.c_str()).what());
    }
    const char * what() const noexcept override { return m_msg.c_str(); }

private:
    std::string m_msg;
};

inline void ThrowException(std::string msg)
{
    GENERIC_THROW(Exception(std::move(msg)))
}

inline void ThrowException(const std::string & msg, int errCode)
{
    GENERIC_THROW(Exception(msg, errCode))
}

}//namespace common
}//namespace generic
#endif//GENERIC_COMMON_EXCEPTION_HPP