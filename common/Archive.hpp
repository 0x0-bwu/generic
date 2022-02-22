/**
 * @file Archive.hpp
 * @author bwu
 * @brief Serialization header files
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_COMMON_ARCHIVE_HPP
#define GENERIC_COMMON_ARCHIVE_HPP
#include "Macros.hpp"

#ifdef BOOST_SERIALIZATION_SUPPORT
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
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/list.hpp>
#include <fstream>
#endif

#endif//GENERIC_COMMON_ARCHIVE_HPP