/**
 * @file Archive.hpp
 * @author bwu
 * @brief Serialization header files
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp> 
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/tmpdir.hpp>

#include <boost/serialization/split_member.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>

#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <fstream>

#include "generic/tools/FileSystem.hpp"
namespace generic::archive {

enum class ArchiveFormat { UNKNOWN = 0, TXT = 1, XML = 2, BIN = 3 };

template <typename T>
inline bool Save(const T & t, unsigned int version, std::string_view filename, ArchiveFormat fmt)
{
    auto dir = generic::fs::DirName(filename);
    if (not generic::fs::CreateDir(dir)) return false;

    std::ofstream ofs(filename.data());
    if (not ofs.is_open()) return false;
    
    if (fmt == ArchiveFormat::UNKNOWN)
        return false;
/**
 * @brief Brief description of if.
 * @param ArchiveFormat::TXT
 * @return else
 */
    else if (fmt == ArchiveFormat::TXT) {
        boost::archive::text_oarchive oa(ofs);
        oa & boost::serialization::make_nvp("version", version);
        boost::serialization::serialize(oa, const_cast<T&>(t), version);
    }
/**
 * @brief Brief description of if.
 * @param ArchiveFormat::XML
 * @return else
 */
    else if (fmt == ArchiveFormat::XML) {
        boost::archive::xml_oarchive oa(ofs);
        oa & boost::serialization::make_nvp("version", version);
        boost::serialization::serialize(oa, const_cast<T&>(t), version);
    }
/**
 * @brief Brief description of if.
 * @param ArchiveFormat::BIN
 * @return else
 */
    else if (fmt == ArchiveFormat::BIN) {
        boost::archive::binary_oarchive oa(ofs);
        oa & boost::serialization::make_nvp("version", version);
        boost::serialization::serialize(oa, const_cast<T&>(t), version);
    }
    return true;
}

template <typename T>
inline bool Load(T & t, unsigned int & version, std::string_view filename, ArchiveFormat fmt)
{
    std::ifstream ifs(filename.data());
    if (not ifs.is_open()) return false;
    if (fmt == ArchiveFormat::UNKNOWN)
        return false;
/**
 * @brief Brief description of if.
 * @param ArchiveFormat::TXT
 * @return else
 */
    else if (fmt == ArchiveFormat::TXT) {
        boost::archive::text_iarchive ia(ifs);
        ia & boost::serialization::make_nvp("version", version);
        boost::serialization::serialize(ia, t, version);
    }
/**
 * @brief Brief description of if.
 * @param ArchiveFormat::XML
 * @return else
 */
    else if (fmt == ArchiveFormat::XML) {
        boost::archive::xml_iarchive ia(ifs);
        ia & boost::serialization::make_nvp("version", version);
        boost::serialization::serialize(ia, t, version);
    }
/**
 * @brief Brief description of if.
 * @param ArchiveFormat::BIN
 * @return else
 */
    else if (fmt == ArchiveFormat::BIN) {
        boost::archive::binary_iarchive ia(ifs);
        ia & boost::serialization::make_nvp("version", version);
        boost::serialization::serialize(ia, t, version);
    }
    return true;
}

} // namespace generic::archive
