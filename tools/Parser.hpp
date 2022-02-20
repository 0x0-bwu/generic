/**
 * @file Parser.hpp
 * @author bwu
 * @brief Parser functions
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_PARSER_HPP
#define GENERIC_PARSER_HPP
#include "generic/common/Exception.hpp"
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/tokenizer.hpp>
#include <sstream>
#include <string>
namespace generic {
///@brief parser functions
namespace parser  {

///@brief parse `input` str by parser `p`
template <typename Parser, typename ... Args>
inline bool Parse(const std::string & input, const Parser & p, Args&& ... args)
{
    auto begin = input.begin(), end = input.end();
    bool ok = boost::spirit::qi::parse(begin, end, p, std::forward<Args>(args)...);
    if(begin != end) return false;
    return ok;
}

///@brief split sting `str` by `seperator`
inline std::vector<std::string> Split(const std::string & str, char seperator)
{
    std::vector<std::string> items;
    using Tokenizer = boost::tokenizer<boost::char_separator<char> >;
    boost::char_separator<char> sep(std::string(1, seperator).c_str());
    Tokenizer tokens(str, sep);
    for(auto iter = tokens.begin(); iter != tokens.end(); ++iter)
        items.push_back(*iter);
    return items;
}

}//namespace parser
}//namespace 
#endif//GENERIC_PARSER_HPP