#ifndef GENERIC_FILESYSTEM_HPP
#define GENERIC_FILESYSTEM_HPP

#if GENERIC_CURRENT_CXX_VERSION >= 20
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
namespace filesystem {

inline std::string CurrentPath();
inline std::string ExecutablePath();
inline bool FileExists(const std::string & filename);
inline bool RemoveFile(const std::string & filename);
inline bool PathExists(const std::string path);
inline bool MakeDir(const std::string & path);
inline bool CreateDir(const std::string & path);
inline std::string DirName(const std::string & path);
inline std::string FileName(const std::string & path);

class FileHelper
{
public:
    explicit FileHelper() = default;

    FileHelper(const FileHelper & ) = delete;
    FileHelper & operator= (const FileHelper & ) = delete;
    ~FileHelper() { Close(); }

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

    bool Reopen(bool truncate, std::string * err = nullptr)
    {
        if(m_filename.empty()){
            if(err) *err = "Fail to reopen file!";
            return false;
        }
        return Open(m_filename, truncate, err);
    }
    
    void Flush()
    {
        std::fflush(m_file);
    }

    void Close()
    {
        if(nullptr != m_file){
            std::fclose(m_file);
            m_file = nullptr;
        }
    }

    void Write(const std::string & buf)
    {
        auto size = buf.size();
        std::fwrite(buf.data(), sizeof(char), size, m_file);
    }

    size_t Size() const
    {
        if(nullptr == m_file) return 0;
        return ::fileno(m_file);
    }

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

// Return true if path exists (file or directory)
inline bool PathExists(const std::string path)
{
    struct stat buffer;
    return (0 == ::stat(path.c_str(), &buffer));
}

// return true on success
inline bool MakeDir(const std::string & path)
{
    return 0 == ::mkdir(path.c_str(), mode_t(0755));
}

// return true on success or if the directory already exists
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

}//namespace filesystem
}//namespace generic
#endif//#define GENERIC_FILESYSTEM_HPP