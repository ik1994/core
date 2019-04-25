#ifndef BIOGEARS_COMMON_FILESYSTEM_PATH_H
#define BIOGEARS_COMMON_FILESYSTEM_PATH_H
/*
    path.h -- A simple class for manipulating paths on Linux/Windows/Mac OS

Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
("Enhancements") to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to the author of this software, without
imposing a separate written license agreement for such Enhancements, then you
hereby grant the following license: a non-exclusive, royalty-free perpetual
license to install, use, modify, prepare derivative works, incorporate into
other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.
*/

//This file has been modified from its original version as part of the biogears project


#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <ShlObj.h>
#include <windows.h>
#else
#include <unistd.h>
#endif
#include <sys/stat.h>

#if defined(__linux)
#include <linux/limits.h>
#endif

namespace biogears {
namespace filesystem {

  /**
 * \brief Simple class for manipulating paths on Linux/Windows/Mac OS
 *
 * This class is just a temporary workaround to avoid the heavy boost
 * dependency until boost::filesystem is integrated into the standard template
 * library at some point in the future.
 */
  class path {
  public:
    using value_type = std::vector<std::string>;
    using iterator = value_type::iterator;
    using const_iterator = value_type::const_iterator;
    using reference = value_type::reference;
    using const_reference = value_type::const_reference;

    enum path_type {
      windows_path = 0,
      posix_path = 1,
#if defined(_WIN32)
      native_path = windows_path
#else
      native_path = posix_path
#endif
    };

    path()
      : m_type(native_path)
      , m_absolute(false)
    {
    }

    path(const path& path)
      : m_type(path.m_type)
      , m_path(path.m_path)
      , m_absolute(path.m_absolute)
    {
    }

    path(path&& path)
      : m_type(path.m_type)
      , m_path(std::move(path.m_path))
      , m_absolute(path.m_absolute)
    {
    }

    path(const char* string) { set(string); }

    path(const std::string& string) { set(string); }

#if defined(_WIN32)
    path(const std::wstring& wstring)
    {
      set(wstring);
    }
    path(const wchar_t* wstring) { set(wstring); }
#endif

    size_t length() const
    {
      return m_path.size();
    }

    auto begin() -> iterator { return m_path.begin(); }

    auto end() -> iterator { return m_path.end(); }

    bool empty() const { return m_path.empty(); }

    bool is_absolute() const { return m_absolute; }

    path make_absolute() const
    {
#if !defined(_WIN32)
      char temp[PATH_MAX];
      if (realpath(str().c_str(), temp) == NULL)
        throw std::runtime_error("Internal error in realpath(): " + std::string(strerror(errno)));
      return path(temp);
#else
      std::wstring value = wstr(), out(MAX_PATH_WINDOWS, '\0');
      DWORD length = GetFullPathNameW(value.c_str(), MAX_PATH_WINDOWS, &out[0], NULL);
      if (length == 0)
        throw std::runtime_error("Internal error in realpath(): " + std::to_string(GetLastError()));
      return path(out.substr(0, length));
#endif
    }

    path make_normal() const
    {
      path working_path;
      for (auto& segment : m_path) {
        if (segment.empty()) {

        } else if (segment == ".") {
          
        } else if (segment == "..")
        {
          working_path = working_path.parent_path();
        } else
        {
          working_path /= segment;
        }
      }
      return working_path;
    }

    bool exists() const
    {
#if defined(_WIN32)
      return GetFileAttributesW(wstr().c_str()) != INVALID_FILE_ATTRIBUTES;
#else
      struct stat sb;
      return stat(str().c_str(), &sb) == 0;
#endif
    }

    size_t file_size() const
    {
#if defined(_WIN32)
      struct _stati64 sb;
      if (_wstati64(wstr().c_str(), &sb) != 0)
        throw std::runtime_error("path::file_size(): cannot stat file \"" + str() + "\"!");
#else
      struct stat sb;
      if (stat(str().c_str(), &sb) != 0)
        throw std::runtime_error("path::file_size(): cannot stat file \"" + str() + "\"!");
#endif
      return (size_t)sb.st_size;
    }

    bool is_directory() const
    {
#if defined(_WIN32)
      DWORD result = GetFileAttributesW(wstr().c_str());
      if (result == INVALID_FILE_ATTRIBUTES)
        return false;
      return (result & FILE_ATTRIBUTE_DIRECTORY) != 0;
#else
      struct stat sb;
      if (stat(str().c_str(), &sb))
        return false;
      return S_ISDIR(sb.st_mode);
#endif
    }

    bool is_file() const
    {
#if defined(_WIN32)
      DWORD attr = GetFileAttributesW(wstr().c_str());
      return (attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY) == 0);
#else
      struct stat sb;
      if (stat(str().c_str(), &sb))
        return false;
      return S_ISREG(sb.st_mode);
#endif
    }

    std::string extension() const
    {
      const std::string& name = filename().string();
      size_t pos = name.find_last_of(".");
      if (pos == std::string::npos)
        return "";
      return name.substr(pos + 1);
    }

    path filename() const
    {
      if (empty())
        return "";
      const std::string& last = m_path[m_path.size() - 1];
      return last;
    }

    path parent_path() const
    {
      path result;
      result.m_absolute = m_absolute;

      if (m_path.empty()) {
        if (!m_absolute)
          result.m_path.push_back("..");
      } else {
        size_t until = m_path.size() - 1;
        for (size_t i = 0; i < until; ++i)
          result.m_path.push_back(m_path[i]);
      }
      return result;
    }

    path& operator/=(const path& other)
    {
      *this = *this / other;
      return *this;
    }
    path operator/(const path& other) const
    {
      if (other.m_absolute)
        throw std::runtime_error("path::operator/(): expected a relative path!");
      if (m_type != other.m_type)
        throw std::runtime_error("path::operator/(): expected a path of the same type!");

      path result(*this);

      for (size_t i = 0; i < other.m_path.size(); ++i)
        result.m_path.push_back(other.m_path[i]);

      return result;
    }

    std::string str(path_type type = native_path) const
    {
      std::ostringstream oss;

      if (m_absolute) {
        if (m_type == posix_path)
          oss << "/";
        else {
          size_t length = 0;
          for (size_t i = 0; i < m_path.size(); ++i)
            // No special case for the last segment to count the NULL character
            length += m_path[i].length() + 1;
          // Windows requires a \\?\ prefix to handle paths longer than MAX_PATH
          // (including their null character). NOTE: relative paths >MAX_PATH are
          // not supported at all in Windows.
          if (length > MAX_PATH_WINDOWS_LEGACY)
            oss << "\\\\?\\";
        }
      }

      for (size_t i = 0; i < m_path.size(); ++i) {
        oss << m_path[i];
        if (i + 1 < m_path.size()) {
          if (type == posix_path)
            oss << '/';
          else
            oss << '\\';
        }
      }

      return oss.str();
    }

    std::string string(path_type type = native_path) const
    {
      return str();
    }

    void set(const std::string& str, path_type type = native_path)
    {
      m_type = type;
      if (type == windows_path) {
        std::string tmp = str;

        // Long windows paths (sometimes) begin with the prefix \\?\. It should only
        // be used when the path is >MAX_PATH characters long, so we remove it
        // for convenience and add it back (if necessary) in str()/wstr().
        static const std::string PREFIX = "\\\\?\\";
        if (tmp.length() >= PREFIX.length()
            && std::mismatch(std::begin(PREFIX), std::end(PREFIX), std::begin(tmp)).first == std::end(PREFIX)) {
          tmp.erase(0, 4);
        }
        m_path = tokenize(tmp, "/\\");
        m_absolute = tmp.size() >= 2 && std::isalpha(tmp[0]) && tmp[1] == ':';
      } else {
        m_path = tokenize(str, "/");
        m_absolute = !str.empty() && str[0] == '/';
      }
    }

    path& operator=(const path& path)
    {
      m_type = path.m_type;
      m_path = path.m_path;
      m_absolute = path.m_absolute;
      return *this;
    }

    path& operator=(path&& path)
    {
      if (this != &path) {
        m_type = path.m_type;
        m_path = std::move(path.m_path);
        m_absolute = path.m_absolute;
      }
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const path& path)
    {
      os << path.str();
      return os;
    }

    bool remove_file()
    {
#if !defined(_WIN32)
      return std::remove(str().c_str()) == 0;
#else
      return DeleteFileW(wstr().c_str()) != 0;
#endif
    }

    bool resize_file(size_t target_length)
    {
#if !defined(_WIN32)
      return ::truncate(str().c_str(), (off_t)target_length) == 0;
#else
      HANDLE handle = CreateFileW(wstr().c_str(), GENERIC_WRITE, 0, nullptr, 0, FILE_ATTRIBUTE_NORMAL, nullptr);
      if (handle == INVALID_HANDLE_VALUE)
        return false;
      LARGE_INTEGER size;
      size.QuadPart = (LONGLONG)target_length;
      if (SetFilePointerEx(handle, size, NULL, FILE_BEGIN) == 0) {
        CloseHandle(handle);
        return false;
      }
      if (SetEndOfFile(handle) == 0) {
        CloseHandle(handle);
        return false;
      }
      CloseHandle(handle);
      return true;
#endif
    }

    static path getcwd()
    {
#if !defined(_WIN32)
      char temp[PATH_MAX];
      if (::getcwd(temp, PATH_MAX) == NULL)
        throw std::runtime_error("Internal error in getcwd(): " + std::string(strerror(errno)));
      return path(temp);
#else
      std::wstring temp(MAX_PATH_WINDOWS, '\0');
      if (!_wgetcwd(&temp[0], MAX_PATH_WINDOWS))
        throw std::runtime_error("Internal error in getcwd(): " + std::to_string(GetLastError()));
      return path(temp.c_str());
#endif
    }

#if defined(_WIN32)
    std::wstring wstr(path_type type = native_path) const
    {
      std::string temp = str(type);
      int size = MultiByteToWideChar(CP_UTF8, 0, &temp[0], (int)temp.size(), NULL, 0);
      std::wstring result(size, 0);
      MultiByteToWideChar(CP_UTF8, 0, &temp[0], (int)temp.size(), &result[0], size);
      return result;
    }

    void set(const std::wstring& wstring, path_type type = native_path)
    {
      std::string string;
      if (!wstring.empty()) {
        int size = WideCharToMultiByte(CP_UTF8, 0, &wstring[0], (int)wstring.size(),
                                       NULL, 0, NULL, NULL);
        string.resize(size, 0);
        WideCharToMultiByte(CP_UTF8, 0, &wstring[0], (int)wstring.size(),
                            &string[0], size, NULL, NULL);
      }
      set(string, type);
    }

    path& operator=(const std::wstring& str)
    {
      set(str);
      return *this;
    }
#endif

    bool operator==(const path& p) const
    {
      return p.m_path == m_path;
    }
    bool operator!=(const path& p) const { return p.m_path != m_path; }

  protected:
    static std::vector<std::string> tokenize(const std::string& string, const std::string& delim)
    {
      std::string::size_type lastPos = 0, pos = string.find_first_of(delim, lastPos);
      std::vector<std::string> tokens;

      while (lastPos != std::string::npos) {
        if (pos != lastPos)
          tokens.push_back(string.substr(lastPos, pos - lastPos));
        lastPos = pos;
        if (lastPos == std::string::npos || lastPos + 1 == string.length())
          break;
        pos = string.find_first_of(delim, ++lastPos);
      }

      return tokens;
    }

  protected:
#if defined(_WIN32)
    static const size_t MAX_PATH_WINDOWS = 32767;
#endif
    static const size_t MAX_PATH_WINDOWS_LEGACY = 260;
    path_type m_type;
    value_type m_path;
    bool m_absolute;
  };

  inline bool create_directory(const path& p)
  {
#if defined(_WIN32)
    return CreateDirectoryW(p.wstr().c_str(), NULL) != 0;
#else
    return mkdir(p.str().c_str(), S_IRWXU) == 0;
#endif
  }

  inline bool create_directories(const path& p)
  {
#if defined(_WIN32)
    return SHCreateDirectory(nullptr, p.make_absolute().wstr().c_str()) == ERROR_SUCCESS;
#else
    if (create_directory(p.str().c_str()))
      return true;

    if (p.empty())
      return false;

    if (errno == ENOENT) {
      if (create_directory(p.parent_path()))
        return mkdir(p.str().c_str(), S_IRWXU) == 0;
      else
        return false;
    }
    return false;
#endif
  }

  inline path absolute(const path p)
  {
    return p.make_absolute();
  }

  inline bool is_directory(const path& p)
  {
    return p.is_directory();
  }
} // namespace filesystem
} //namespace biogears

#endif //BIOGEARS_COMMON_FILESYSTEM_PATH_H
