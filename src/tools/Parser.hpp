/**
 * @file Parser.hpp
 * @author bwu
 * @brief Parser functions
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "generic/common/Exception.hpp"
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
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

///@brief print line number and line content if error is met 
template <typename Iterator>
struct ErrorHandler
{
	template <typename, typename, typename>
	struct result { using type = void; };

	ErrorHandler(Iterator first, Iterator last, std::ostream & os = std::cout)
	  : first(first), last(last), os(os) {}

	template <typename Message, typename What>
	void operator() (const Message & message, const What & what, Iterator errPos) const
	{
		int line;
		Iterator lineStart = GetPos(errPos, line);
		if (errPos != last) {
			os << message << what << " line " << line << ':' << std::endl;
			os << GetLine(lineStart) << std::endl;
			for (; lineStart != errPos; ++lineStart)
				os << ' ';
			os << '^' << std::endl;
		}
		else {
			os << "Unexpected end of content. ";
			os << message << what << " line " << line << std::endl;
		}
	}

	Iterator GetPos(Iterator errPos, int & line) const
	{
		line = 1;
		Iterator i = first;
		Iterator lineStart = first;
		while (i != errPos) {
			bool eol = false;
			if (i != errPos && *i == '\r') {// CR
				eol = true;
				lineStart = ++i;
			}
			if (i != errPos && *i == '\n') {// LF
				eol = true;
				lineStart = ++i;
			}
			if (eol) ++line;
			else ++i;
		}
		return lineStart;
	}

	std::string GetLine(Iterator errPos) const
	{
		Iterator i = errPos;
		// position i to the next EOL
		while (i != last && (*i != '\r' && *i != '\n')) { ++i; }
		return std::string(errPos, i);
	}

	Iterator first;
	Iterator last;
	std::vector<Iterator> iters;
    std::ostream & os;
};

} // namespace parser
} // namespace generic