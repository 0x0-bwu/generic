#ifndef GENERIC_PROGRAMOPTIONS_HPP
#define GENERIC_PROGRAMOPTIONS_HPP
#include "generic/common/Exception.hpp"
#include <type_traits>
#include <utility>
#include <sstream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
namespace generic {
namespace program_options {

enum class Argument { No = 0, Required, Optional };
enum class OptionName { Unspecified, ShortName, LongName };
enum class Attribute { Inactive = 0, Hidden = 1, Required = 2, Optional = 3, Advanced = 4, Expert = 5 };

class Option
{
    friend class OptionParser;
public:
    Option(std::string shortName, std::string longName, std::string description)
     : m_shortName(std::move(shortName)), m_longName(std::move(longName)), m_description(std::move(description)) {}
    virtual ~Option() = default;

    char ShortName() const
    {
        if(!m_shortName.empty())
            return m_shortName.front();
        return 0;
    }

    const std::string & LongName() const
    {
        return m_longName;
    }

    std::string Name(OptionName whatName, bool withHypen = false) const
    {
        if(whatName == OptionName::ShortName)
            return m_shortName.empty() ? "" : ((withHypen ? "-" : "") + m_shortName);
        if(whatName == OptionName::LongName)
            return m_longName.empty()  ? "" : ((withHypen ? "--" : "") + m_longName);
        return "";
    }

    const std::string & Description() const
    {
        return m_description;
    }

    void SetAttribute(Attribute attribute)
    {
        m_attribute = attribute;
    }

    Attribute GetAttribute() const
    {
        return m_attribute;
    }
    
    virtual bool GetDefault(std::ostream & out) const = 0;
    virtual Argument ArgumentType() const = 0;
    virtual size_t Count() const = 0;
    virtual bool isSet() const = 0;

protected:
    virtual void Parse(OptionName whatName, const char * value) = 0;
    virtual void Clear() = 0;

protected:
    std::string m_shortName;
    std::string m_longName;
    std::string m_description;
    Attribute m_attribute;
};

struct InvalidOption : public std::invalid_argument
{
    enum class Error
    {
        InvalidArgument,
        MissingArgument,
        TooManyArgument,
        MissingOption
    };

    const Option * option;
    Error error;
    OptionName whatName;
    std::string value;

    InvalidOption(const Option * _option, InvalidOption::Error _error, OptionName _whatName, std::string _value, const std::string & text)
     : std::invalid_argument(text.c_str()), option(_option), error(_error), whatName(_whatName), value(std::move(_value)) {}

    InvalidOption(const Option * option, InvalidOption::Error error, const std::string & text)
     : InvalidOption(option, error, OptionName::Unspecified, "", text) {}
};

//Value option with optional default value
template <typename T>
class Value : public Option
{
public:
    Value(std::string shortName, std::string longName, std::string description)
     : Option(std::forward<std::string>(shortName),
              std::forward<std::string>(longName),
              std::forward<std::string>(description)) {}
    
    Value(std::string shortName, std::string longName, std::string description, const T & defaultValue, T * assignTo = nullptr)
     : Option(std::forward<std::string>(shortName),
              std::forward<std::string>(longName),
              std::forward<std::string>(description)) { m_assignTo = assignTo; SetDefault(defaultValue); }

    virtual ~Value() = default;

    size_t Count() const override
    {
        return m_values.size();
    }

    bool isSet() const override
    {
        return !m_values.empty();
    }

    void AssignTo(T * var)
    {
        m_assignTo = var;
        UpdateReference();
    }

    void SetValue(const T & value)
    {
        Clear();
        AddValue(value);
    }

    T GetValue(size_t index = 0) const
    {
        if(!isSet() && m_default) return *m_default;

        if(!isSet() || index >= Count()){

            std::stringstream err;
            if(!isSet()) err << "option no set: \"";
            else err << "index out of range (" << index << ") for \"";

            if(ShortName() != 0) err << "-" << ShortName();
            else err << "--" << LongName();

            err << "\"";
            GENERIC_THROW(std::out_of_range(err.str()))
        }

        return m_values[index];
    }

    T ValueOr(const T & defaultValue, size_t index) const
    {
        if(index < m_values.size()) return m_values[index];
        else if(m_default) return *m_default;
        else return defaultValue;
    }

    void SetDefault(const T & value)
    {
        m_default.reset(new T(value));
        UpdateReference();
    }

    bool hasDefault() const
    {
        return m_default != nullptr;
    }

    T GetDefault() const
    {
        if(!hasDefault())
            GENERIC_THROW(std::runtime_error("no default value set"))
        return m_default;
    }

    bool GetDefault(std::ostream & out) const override
    {
        if(!hasDefault()) return false;
        out << *m_default;
        return true;
    }

    Argument ArgumentType() const override
    {
        return Argument::Required;
    }

protected:
    void Parse(OptionName whatName, const char * value) override
    {
        T parsedVal;
        std::string strVal;
        if(value != nullptr)
            strVal = value;
        
        std::istringstream is(strVal);
        int valRead = 0;
        while(is.good()){
            if(is.peek() != EOF)
                is >> parsedVal;
            else break;
            valRead++;
        }

        if(is.fail()){
            std::string err = "invalid argument for " + Name(whatName, true) + ": '" + strVal + "'";
            GENERIC_THROW(InvalidOption(this, InvalidOption::Error::InvalidArgument, whatName, value, err))
        }

        if(valRead > 1){
            std::string err = "too many arguments for " + Name(whatName, true) + ": '" + strVal + "'";
            GENERIC_THROW(InvalidOption(this, InvalidOption::Error::TooManyArgument, whatName, value, err))
        }

        if(strVal.empty()){
            std::string err = "missing argument for " + Name(whatName, true);
            GENERIC_THROW(InvalidOption(this, InvalidOption::Error::MissingArgument, whatName, "", err))
        }

        this->AddValue(parsedVal);
    }

    virtual void UpdateReference()
    {
        if(m_assignTo){
            if(!isSet() && m_default)
                *m_assignTo = *m_default;
            else if(isSet())
                *m_assignTo = m_values.back();
        }
    }

    virtual void AddValue(const T & value)
    {
        m_values.push_back(value);
        UpdateReference();
    }

    void Clear() override
    {
        m_values.clear();
        UpdateReference();
    }

protected:
    std::unique_ptr<T> m_default;
    T * m_assignTo = nullptr;
    std::vector<T> m_values;
};

//Value option with implicit default value
template <typename T>
class Implicit : public Value<T>
{
public:
    Implicit(std::string shortName, std::string longName, std::string description, const T & implicitVal, T * assignTo = nullptr)
     : Value<T>(std::forward<std::string>(shortName),
                std::forward<std::string>(longName),
                std::forward<std::string>(description), implicitVal, assignTo) {}
    
    Argument ArgumentType() const override { return Argument::Optional; }

protected:
    void Parse(OptionName whatName, const char * value) override
    {
        if((value != nullptr) && (std::char_traits<char>::length(value) > 0))
            Value<T>::Parse(whatName, value);
        else Value<T>::AddValue(*Value<T>::m_default);
    }
};

//Value option without value
class Switch : public Value<bool>
{
public:
    Switch(std::string shortName, std::string longName, std::string description, bool * assignTo = nullptr)
     : Value<bool>(std::forward<std::string>(shortName),
                   std::forward<std::string>(longName),
                   std::forward<std::string>(description), false, assignTo){}
    
    void SetDefault(const bool & value) = delete;
    Argument ArgumentType() const override { return Argument::No; }

protected:
    void Parse(OptionName , const char * ) override
    {
        AddValue(true);
    }
};

using OptionPtr = std::shared_ptr<Option>;
class OptionParser;
class OptionPrinter
{
public:
    explicit OptionPrinter(const OptionParser * parser) : m_parser(parser) {}
    virtual ~OptionPrinter() = default;

    virtual std::string Print(const Attribute & showLevel = Attribute::Optional) const = 0;

protected:
    const OptionParser * m_parser;
};

class ConsoleOptionPrinter : public OptionPrinter
{
public:
    explicit ConsoleOptionPrinter(const OptionParser * parser) : OptionPrinter(parser) {}
    virtual ~ConsoleOptionPrinter() = default;

    std::string Print(const Attribute & showLevel = Attribute::Optional) const override;
private:
    std::string toString(OptionPtr option) const;
};


class OptionParser
{
public:
    explicit OptionParser(std::string description = "") : m_description(std::move(description)) {}
    virtual ~OptionParser() = default;

    template <typename T, typename... Args>
    std::shared_ptr<T> Add(Args&&... args)
    {
        return Add<T, Attribute::Optional>(std::forward<Args>(args)...);
    }

    template <typename T, Attribute attribute, typename... Args>
    std::shared_ptr<T> Add(Args&&... args)
    {
        static_assert(std::is_base_of<Option, typename std::decay<T>::type>::value);
        
        auto option = std::make_shared<T>(std::forward<Args>(args)...);

        for(const auto & op : m_options){
            if((option->ShortName() != 0) && (option->ShortName() == op->ShortName())){
                std::string err = "duplicate short option name '-" + std::string(1, option->ShortName()) + "'";
                GENERIC_THROW(std::invalid_argument(err))
            }
            if(!option->LongName().empty() && (option->LongName() == op->LongName())){
                std::string err = "duplicate long option name '--" + option->LongName() + "'";
                GENERIC_THROW(std::invalid_argument(err))
            }
        }
        option->SetAttribute(attribute);
        m_options.push_back(option);
        return option;
    }

    bool Parse(const std::string & configFile)
    {
        auto trim = [](std::string & s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) { return !std::isspace(ch); }));
            s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch){ return !std::isspace(ch); }).base(), s.end());
            return s;
        };

        auto trimCopy = [trim](const std::string & s) { std::string copy(s); return trim(copy); };

        auto split = [trimCopy](const std::string & s) -> std::pair<std::string, std::string> {
            auto pos = s.find('=');
            if(pos == std::string::npos) return {"",""};
            return {trimCopy(s.substr(0, pos)), trimCopy(s.substr(pos + 1, std::string::npos))};
        };

        std::ifstream file(configFile.c_str());
        if(!file.is_open()) return false;

        std::string line, section;
        while(std::getline(file, line)){
            trim(line);
            if(line.empty()) continue;
            if(line.front() == '#') continue;
            if(line.front() == '[' && line.back() == ']'){
                section = trimCopy(line.substr(1, line.size() - 2));
                continue;
            }

            auto keyVal = split(line);
            if(keyVal.first.empty()) continue;
            auto key = section.empty() ? keyVal.first : section + "." + keyVal.first;
            auto option = FindOption(key);
            if(option && option->GetAttribute() == Attribute::Inactive) option = nullptr;
            if(option) option->Parse(OptionName::LongName, keyVal.second.c_str());
            else m_unknowOptions.push_back(key);
        }
        return true;
    }

    void Parse(int argc, const char * const argv[])
    {
        for(auto n = 1; n < argc; ++n){
            const std::string arg(argv[n]);
            if(arg == "--"){
                for(int m = n + 1; m < argc; ++m)
                    m_noOptionArgs.emplace_back(argv[m]);
            }
            else if(arg.find("--") == 0){
                std::string opt = arg.substr(2);
                std::string optArg;
                auto equalIdx = opt.find('=');
                if(equalIdx != std::string::npos){
                    optArg = opt.substr(equalIdx + 1);
                    opt.resize(equalIdx);
                }

                auto option = FindOption(opt);
                if(option && option->GetAttribute() == Attribute::Inactive) option = nullptr;
                if(option){
                    if(option->ArgumentType() == Argument::No){
                        if(!optArg.empty()) option = nullptr;
                    }
                    else if(option->ArgumentType() == Argument::Required){
                        if(optArg.empty() && n < argc - 1) optArg = argv[++n];
                    }
                }

                if(option) option->Parse(OptionName::LongName, optArg.c_str());
                else m_unknowOptions.push_back(arg);
            }
            else if(arg.find('-') == 0){
                std::string opt = arg.substr(1);
                bool unknow = false;
                for(size_t m = 0; m < opt.size(); ++m){
                    char c = opt[m];
                    std::string optArg;
                    auto option = FindOption(c);
                    if(option && option->GetAttribute() == Attribute::Inactive) option = nullptr;
                    if(option){
                        if(option->ArgumentType() == Argument::Required){
                            optArg = opt.substr(m + 1);
                            if(optArg.empty() && n < argc - 1) optArg = argv[++n];
                            m = opt.size();
                        }
                        else if(option->ArgumentType() == Argument::Optional){
                            optArg = opt.substr(m + 1);
                            m = opt.size();
                        }
                    }
                    if(option) option->Parse(OptionName::ShortName, optArg.c_str());
                    else unknow = true;
                }
                if(unknow) m_unknowOptions.push_back(arg);
            }
            else m_noOptionArgs.push_back(arg);
        }

        for(auto & opt : m_options){
            if(opt->GetAttribute() == Attribute::Required && !opt->isSet()){
                std::string option = opt->LongName().empty() ? std::string(1, opt->ShortName()) : opt->LongName();
                std::string error = "Error: option \"" + option + "\" is required";
                GENERIC_THROW(InvalidOption(opt.get(), InvalidOption::Error::MissingOption, error))
            }
        }
    }

    void Reset()
    {
        m_unknowOptions.clear();
        m_noOptionArgs.clear();
        for(auto & opt : m_options)
            opt->Clear();
    }

    std::string Help(Attribute showLevel = Attribute::Optional) const
    {
        ConsoleOptionPrinter printer(this);
        return printer.Print(showLevel);
    }

    std::string Description() const { return m_description; }

    const std::vector<OptionPtr> & Options() const { return m_options; }

    const std::vector<std::string> & NonOptionArgs() const { return m_noOptionArgs; }

    const std::vector<std::string> & UnknowOptions() const { return m_unknowOptions; }

    template <typename T>
    std::shared_ptr<T> GetOption(const std::string & longName) const
    {
        auto option = FindOption(longName);
        return TryGetOption<T>(option, longName);
    }

    template <typename T>
    std::shared_ptr<T> GetOption(char shortName) const
    {
        auto option = FindOption(shortName);
        return TryGetOption<T>(option, std::string(1, shortName));
    }

protected:
    OptionPtr FindOption(const std::string & longName) const
    {
        for(const auto & option : m_options){
            if(option->LongName() == longName) return option;
        }
        return nullptr;
    }

    OptionPtr FindOption(char shortName) const
    {
        for(const auto & option : m_options){
            if(option->ShortName() == shortName) return option;
        }
        return nullptr;
    }

    template <typename T>
    std::shared_ptr<T> TryGetOption(OptionPtr option, const std::string & name) const
    {
        if(!option){
            std::string err = "option not found: " + name;
            GENERIC_THROW(std::invalid_argument(err))
        }
        auto result = std::dynamic_pointer_cast<T>(option);
        if(!result){
            std::string err = "cannot cast option to T: " + name;
            GENERIC_THROW(std::invalid_argument(err))
        }
        return result;
    }
    
protected:
    std::vector<OptionPtr> m_options;
    std::string m_description;
    std::vector<std::string> m_noOptionArgs;
    std::vector<std::string> m_unknowOptions;
};

inline std::string ConsoleOptionPrinter::Print(const Attribute & showLevel) const
{
    if(nullptr == m_parser) return std::string{};
    if(showLevel < Attribute::Optional) return std::string{};

    std::stringstream ss;
    if(!m_parser->Description().empty()){
        ss << m_parser->Description() << ':' << GENERIC_DEFAULT_EOL;
    }

    size_t optRightMargin(20);
    size_t maxDescriptionLeftMargin(40);

    for(const auto & option : m_parser->Options())
        optRightMargin = std::max(optRightMargin, toString(option).size() + 2);
    optRightMargin = std::min(maxDescriptionLeftMargin - 2, optRightMargin);

    for(const auto & option : m_parser->Options()){
        if(option->GetAttribute() <= Attribute::Hidden ||
            option->GetAttribute() > showLevel) continue;
        
        auto optStr = toString(option);
        if(optStr.size() < optRightMargin)
            optStr.resize(optRightMargin, ' ');
        else{
            optStr += "\n" + std::string(optRightMargin, ' ');
        }
        ss << optStr;

        std::string line;
        std::vector<std::string> lines;
        std::stringstream description(option->Description());
        while(std::getline(description, line, '\n'))
            lines.push_back(line);
        
        std::string empty(optRightMargin, ' ');
        for(size_t i = 0; i < lines.size(); ++i){
            if(i > 0) ss << GENERIC_DEFAULT_EOL << empty;
            ss << lines.at(i);
        }
        ss << std::endl;
    }

    return ss.str();
}

inline std::string ConsoleOptionPrinter::toString(OptionPtr option) const
{
    std::stringstream ss;
    if(option->ShortName() != 0){
        ss << "  -" << option->ShortName();
        if(!option->LongName().empty())
            ss << ", ";
    }
    else ss << "  ";
    if(!option->LongName().empty())
        ss << "--" << option->LongName();
    
    if(option->ArgumentType() == Argument::Required){
        ss << " arg";
        std::stringstream defaultStr;
        if(option->GetDefault(defaultStr)){
            if(!defaultStr.str().empty()){
                ss << " (="<< defaultStr.str() << ")";
            }
        }
    }
    else if(option->ArgumentType() == Argument::Optional){
        std::stringstream defaultStr;
        if(option->GetDefault(defaultStr)){
            ss <<" [=arg(=" << defaultStr.str() << ")]";
        }
    }
    return ss.str();
}

}//namespace program_options
}//namespace generic

namespace {
inline std::ostream & operator<< (std::ostream & os, const generic::program_options::OptionParser & op)
{
    return os << op.Help();
}
}
#endif//GENERIC_PROGRAMOPTIONS_HPP