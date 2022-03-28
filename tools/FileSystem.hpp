/**
 * @file FileSystem.hpp
 * @author bwu
 * @brief File system related functions
 * @version 0.1
 * @date 2022-02-22
 */
#ifndef GENERIC_FILESYSTEM_HPP
#define GENERIC_FILESYSTEM_HPP

#if GENERIC_CURRENT_CXX_VERSION >= 17
#include <filesystem>
#endif

#ifdef GENERIC_OS_WINDOWS
#include <windows.h>
#else
#include <limits.h>
#include <unistd.h>
#endif

#include "Tools.hpp"
#include <fstream>
#include <cstdio>
namespace generic{
///@brief filesystem related functions
namespace filesystem {

///@brief returns absolute path of the current working directory
inline std::string CurrentPath();

///@brief returns absolute path of the current executable binary
inline std::string ExecutablePath();

///@brief checks if file exist
inline bool FileExists(const std::string & filename);

///@brief removes a file from disk
inline bool RemoveFile(const std::string & filename);

///@brief checks if path exists (file or directory)
inline bool PathExists(const std::string path);

///@brief makes a directory of given path
inline bool MakeDir(const std::string & path);

///@brief makes a directory of given path, return true on success or if the directory already exists
inline bool CreateDir(const std::string & path);

///@brief gets folder name of given file path
inline std::string DirName(const std::string & path);

///@brief gets file name of given file path
inline std::string FileName(const std::string & path);

///@brief gets base file name of given file path
inline std::string BaseName(const std::string & path);

#ifndef GENERIC_OS_WINDOWS
///@brief checks if given folder is writable
inline bool isDirWritable(const std::string & path);
#endif//GENERIC_OS_WINDOWS

#if GENERIC_CURRENT_CXX_VERSION >= 17
///@brief empties given folder
inline bool RemoveDir(const std::string & path);
#endif

///@brief represents a helper class for file writing
class FileHelper
{
public:
    explicit FileHelper() = default;

    FileHelper(const FileHelper & ) = delete;
    FileHelper & operator= (const FileHelper & ) = delete;
    ~FileHelper() { Close(); }

    /**
     * @brief open a file on disk for write
     * @param[in] filename file to open
     * @param[in] truncate whether truncate when openning
     * @param[out] err error message if filed to open file
     * @return whether the file opened
     */
    bool Open(const std::string & filename, bool truncate = true, std::string * err = nullptr)
    {
        Close();
        m_filename = filename;

        auto * mode = "ab";
        auto * truncMode = "wb";
        for(auto t = 0; t < m_openTries; ++t){
            CreateDir(DirName(filename));
            if(truncate){
                std::FILE * tmp = std::fopen(filename.c_str(), truncMode);
                if(!tmp) continue;
                std::fclose(tmp);
            }

            m_file = std::fopen(filename.c_str(), mode);
            if(m_file) return true;

            tools::SleepMilliseconds(m_openInterval);
        }
        if(err) *err = "Fail to open file " + filename + " for writing!";
        return false;
    }

    ///@brief reopens the closed file
    bool Reopen(bool truncate, std::string * err = nullptr)
    {
        if(m_filename.empty()){
            if(err) *err = "Fail to reopen file!";
            return false;
        }
        return Open(m_filename, truncate, err);
    }
    
    ///@brief flushes the string stream in cache to file
    void Flush()
    {
        std::fflush(m_file);
    }

    ///@brief closes current file
    void Close()
    {
        if(nullptr != m_file){
            std::fclose(m_file);
            m_file = nullptr;
        }
    }

    ///@brief writes the buf to cache
    void Write(const std::string & buf)
    {
        auto size = buf.size();
        std::fwrite(buf.data(), sizeof(char), size, m_file);
    }

    ///@brief returns current file size
    size_t Size() const
    {
        if(nullptr == m_file) return 0;
        return ::fileno(m_file);
    }

    ///@brief gets current file name
    const std::string & Filename() const
    {
        return m_filename;
    }

private:
    const int m_openTries = 5;
    const unsigned int m_openInterval = 10;
    std::FILE * m_file{nullptr};
    std::string m_filename;
};

inline std::string CurrentPath()
{
#if GENERIC_CURRENT_CXX_VERSION >= 20
    return std::filesystem::current_path();
#else
    char buffer[1024];
    char *cwd = getcwd(buffer, sizeof(buffer));
    return std::string(cwd);
#endif
}

inline std::string ExecutablePath()
{
#ifdef GENERIC_OS_WINDOWS
    wchar_t path[4096] = { 0 };
    GetModuleFileNameW(NULL, path, MAX_PATH);
    std::wstring ws(path);
    return std::string(ws.begin(), ws.end());
#else
    char result[4096];
    auto count = readlink("/proc/self/exe", result, 4096);
    return std::string(result, (count > 0) ? count : 0);
#endif
}

inline bool FileExists(const std::string & filename)
{
    std::ifstream f(filename);
    return f.good();
}

inline bool RemoveFile(const std::string & filename)
{
    return 0 == std::remove(filename.c_str());
}

inline bool PathExists(const std::string path)
{
    struct stat buffer;
    return (0 == ::stat(path.c_str(), &buffer));
}

inline bool MakeDir(const std::string & path)
{
    return 0 == ::mkdir(path.c_str(), mode_t(0755));
}

inline bool CreateDir(const std::string & path)
{
    if(PathExists(path)) return true;
    if(path.empty()) return false;
    size_t offset = 0;
    do {
        auto tokenPos = path.find_first_of(GENERIC_FOLDER_SEPS, offset);
        if(tokenPos == std::string::npos) tokenPos = path.size();
        auto subdir = path.substr(0, tokenPos);
        if(!subdir.empty() && !PathExists(subdir) && !MakeDir(subdir)) return false;
        offset = tokenPos + 1;
    } while (offset < path.size());

    return true;
}

inline std::string DirName(const std::string & path)
{
    auto pos = path.find_last_of(GENERIC_FOLDER_SEPS);
    return pos != std::string::npos ? path.substr(0, pos) : std::string{};
}

inline std::string FileName(const std::string & path)
{
    auto pos = path.find_last_of(GENERIC_FOLDER_SEPS);
    return pos != std::string::npos ? path.substr(pos + 1) : path.substr(0);
}

inline std::string BaseName(const std::string & path){
    auto filename = FileName(path);
    auto pos = filename.find_first_of('.');
    return pos != std::string::npos ? filename.substr(0, pos) : filename.substr(0);
}

#ifndef GENERIC_OS_WINDOWS
inline bool isDirWritable(const std::string & path)
{
    if(access(path.c_str(), W_OK) == 0) return true;
    return false;
}
#endif//GENERIC_OS_WINDOWS

#if GENERIC_CURRENT_CXX_VERSION >= 17
inline bool RemoveDir(const std::string & path)
{
    return std::filesystem::remove_all(path) > 0;
}
#endif

}//namespace filesystem
}//namespace generic
#endif//GENERIC_FILESYSTEM_HPP