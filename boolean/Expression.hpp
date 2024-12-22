/**
 * @file Expression.hpp
 * @author bwu
 * @brief boolean expresion
 * @version 0.1
 * @date 2024-12-21
 */
#pragma once
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/variant/recursive_wrapper.hpp>
#include <boost/spirit/include/phoenix.hpp>
#include <boost/spirit/include/qi.hpp>
#include <string_view>
namespace generic::boolean::expr {

namespace phx = boost::phoenix;
namespace qi = boost::spirit::qi;

struct Or;
struct And;
struct Not;
struct Xor;

using Variable = std::string;

template <typename Tag> struct UnaryOperation;
template <typename Tag> struct BinaryOperation;

using OperationOr = BinaryOperation<Or>;
using OperationNot = UnaryOperation<Not>;
using OperationAnd = BinaryOperation<And>;
using OperationXor = BinaryOperation<Xor>;

using Expression = boost::variant<Variable,
                   boost::recursive_wrapper<OperationOr>,
                   boost::recursive_wrapper<OperationNot>,
                   boost::recursive_wrapper<OperationAnd>,
                   boost::recursive_wrapper<OperationXor>>;

template <typename Tag>
struct UnaryOperation
{
    Expression expr;
    explicit UnaryOperation(Expression e)
     : expr(std::move(e)) {}
};

template <typename Tag>
struct BinaryOperation
{
    Expression left, right;
    explicit BinaryOperation(Expression l, Expression r)
     : left(std::move(l)), right(std::move(r)) {}
};

struct ExpressionPrinter : boost::static_visitor<void>
{
    explicit ExpressionPrinter(std::ostream & os) : m_os(os) {}
    void operator()(const Variable & v) const { m_os << v; }
    void operator()(const OperationOr  & op) const { Print(" | ", op.left, op.right); }
    void operator()(const OperationAnd & op) const { Print(" & ", op.left, op.right); }
    void operator()(const OperationXor & op) const { Print(" ^ ", op.left, op.right); }
    void operator()(const OperationNot & op) const
    {
        m_os << "(!";
        boost::apply_visitor(*this, op.expr);
        m_os << ')';
    }
private:
    void Print(const std::string & op, const Expression & l, const Expression & r) const
    {
        m_os << '(';
        boost::apply_visitor(*this, l);
        m_os << op;
        boost::apply_visitor(*this, r);
        m_os << ')';
    }

private:
    std::ostream & m_os;
};

template <typename Iterator, typename Skipper = qi::space_type>
struct ExpressionGrammar : qi::grammar<Iterator, Expression(), Skipper> 
{
    qi::symbols<char, char> escaped;
    qi::rule<Iterator, Variable(), Skipper> variable;
    qi::rule<Iterator, Expression(), Skipper> opAnd, opNot, opXor, opOr, simple, expression;
    ExpressionGrammar() : ExpressionGrammar::base_type(expression)
    {
        using qi::_1;
        using qi::lit;
        using qi::_val;
        using qi::char_;
        escaped.add("\\/", '/')("\\[", '[')("\\]", ']');
        expression = opOr.alias();
        opNot = (lit('!') > simple)[_val = phx::construct<OperationNot>(_1)] | simple[_val = _1];
        opOr  = opXor[_val = _1] >> *(lit('|') >> opXor[_val = phx::construct<OperationOr >(_val, _1)]);
        opXor = opAnd[_val = _1] >> *(lit('^') >> opAnd[_val = phx::construct<OperationXor>(_val, _1)]);
        opAnd = opNot[_val = _1] >> *(lit('&') >> opNot[_val = phx::construct<OperationAnd>(_val, _1)]);
        simple = ('(' > expression > ')') | variable;
        variable = char_("A-Za-z_") >> *(escaped | (char_("A-Za-z-0-9_[]") - '\\'));

        BOOST_SPIRIT_DEBUG_NODE(opOr);
        BOOST_SPIRIT_DEBUG_NODE(opXor);
        BOOST_SPIRIT_DEBUG_NODE(opAnd);
        BOOST_SPIRIT_DEBUG_NODE(opNot);
        BOOST_SPIRIT_DEBUG_NODE(simple);
        BOOST_SPIRIT_DEBUG_NODE(variable);
        BOOST_SPIRIT_DEBUG_NODE(expression);
    }
};

template <typename Iterator>
inline bool ParseExpression(Iterator begin, Iterator end, Expression & expr)
{
    ExpressionGrammar<Iterator> grammar;
    auto res = qi::phrase_parse(begin, end, grammar, qi::space, expr);
    return res and (begin == end);
}

inline bool ParseExpression(std::string_view sv, Expression & expr)
{
    return ParseExpression(sv.begin(), sv.end(), expr);
}

} // namespace generic::boolean::expr

namespace {

inline std::ostream & operator<< (std::ostream & os, const generic::boolean::expr::Expression & e)
{
    boost::apply_visitor(generic::boolean::expr::ExpressionPrinter(os), e);
    return os;
}

}