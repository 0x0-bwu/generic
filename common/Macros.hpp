/**
 * @file Macros.hpp
 * @author bwu
 * @brief Macro defines
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once
#include "Version.hpp"

#define BOOST_SERIALIZATION_SUPPORT 1
#define BOOST_GIL_IO_PNG_SUPPORT    1
#define EIGEN_LIBRARY_SUPPORT       1
#define GENERIC_NO_EXCEPTION        0    
#define GENERIC_NO_ASSERT           0

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
    #define GENERIC_OS_WINDOWS
        #ifdef __WIN64
            #define GENERIC_OS_WINDOWS_64
        #else
            #define GNEERIC_OS_WINDOWS_32
        #endif
#elif defined(__APPLE__)
    #include <TargetConditionals.h>
    #if TARGET_OS_MAC
        #define GENERIC_OS_MAC
    #endif
#else
    #define GENERIC_OS_LINUX
#endif

#ifdef GENERIC_OS_WINDOWS
    #define GENERIC_FOLDER_SEPS "\\/"
    #define GENERIC_DEFAULT_EOL "\r\n"
#else
    #define GENERIC_FOLDER_SEPS "/"
    #define GENERIC_DEFAULT_EOL "\n"
#endif
