/**
 * @file FileSystem.hpp
 * @author bwu
 * @brief File system related functions
 * @version 0.1
 * @date 2022-02-22
 */
#pragma once

#include "Tools.hpp"
#include <filesystem>
namespace generic{
///@brief filesystem related functions
namespace fs {

///@brief returns absolute path of the current working directory
inline std::filesystem::path CurrentPath();

///@brief returns absolute path of the current executable binary
inline std::filesystem::path ExecutablePath();

///@brief checks if file exist
inline bool FileExists(const std::filesystem::path & filename);

///@brief removes a file from disk
inline bool RemoveFile(const std::filesystem::path & filename);

///@brief checks if path exists (file or directory)
inline bool PathExists(const std::filesystem::path & path);

///@brief makes a directory of given path, return true on success or if the directory already exists
inline bool CreateDir(const std::filesystem::path & path);

///@brief gets folder name of given file path
inline std::filesystem::path DirName(const std::filesystem::path & path);

///@brief gets file name of given file path
inline std::filesystem::path FileName(const std::filesystem::path & path);

///@brief gets base file name of given file path
inline std::filesystem::path BaseName(const std::filesystem::path & path);


///@brief empties given folder
inline bool RemoveDir(const std::filesystem::path & path);

///@brief gets parent folder of given file path
inline std::filesystem::path ParentPath(const std::filesystem::path & path);

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
    bool Open(std::string_view filename, bool truncate = true, std::string * err = nullptr)
    {
        Close();
        m_filename = filename;

        auto * mode = "ab";
        auto * truncMode = "wb";
        for(auto t = 0; t < m_openTries; ++t){
            CreateDir(DirName(filename));
            if(truncate){
                std::FILE * tmp = std::fopen(filename.data(), truncMode);
                if(!tmp) continue;
                std::fclose(tmp);
            }

            m_file = std::fopen(filename.data(), mode);
            if(m_file) return true;

            tools::SleepMilliseconds(m_openInterval);
        }
        if(err) *err = "Fail to open file " + std::string(filename) + " for writing!";
        return false;
    }

    ///@brief reopens the closed file
    bool Reopen(bool truncate, std::string * err = nullptr)
    {
        if (m_filename.empty()){
            if (err) *err = "Fail to reopen file!";
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
    void Write(std::string_view buf)
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

inline std::filesystem::path CurrentPath()
{
#if GENERIC_CURRENT_CXX_VERSION >= 20
    return std::filesystem::current_path();
#else
    char buffer[1024];
    char *cwd = getcwd(buffer, sizeof(buffer));
    return std::filesystem::path(cwd);
#endif
}

inline std::filesystem::path ExecutablePath()
{
#ifdef GENERIC_OS_WINDOWS
    wchar_t path[4096] = { 0 };
    GetModuleFileNameW(NULL, path, MAX_PATH);
    std::wstring ws(path);
    return std::filesystem::path(ws.begin(), ws.end());
#else
    char result[4096];
    auto count = readlink("/proc/self/exe", result, 4096);
    auto path = std::string(result, (count > 0) ? count : 0);
    return std::filesystem::path(path);
#endif
}

inline bool FileExists(const std::filesystem::path & filename)
{
    return std::filesystem::exists(filename);
}

inline bool RemoveFile(const std::filesystem::path & filename)
{
    return std::filesystem::remove(filename);
}

inline bool PathExists(const std::filesystem::path & path)
{
    return std::filesystem::exists(path);
}

inline bool CreateDir(const std::filesystem::path & path)
{
    if (PathExists(path)) return true;
    return std::filesystem::create_directories(path);
}

inline std::filesystem::path DirName(const std::filesystem::path & path)
{
    return path.parent_path();
}

inline std::filesystem::path FileName(const std::filesystem::path & path)
{
    return path.filename();
}

inline std::filesystem::path BaseName(const std::filesystem::path & path)
{
    auto p = path;
    while (not p.extension().empty()) p = p.stem();
    return p;
}

inline bool RemoveDir(const std::filesystem::path & path)
{
    return std::filesystem::remove_all(path) > 0;
}

inline std::filesystem::path ParentPath(const std::filesystem::path & path)
{
    return path.parent_path();
}

}//namespace fs
}//namespace generic