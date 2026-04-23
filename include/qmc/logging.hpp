// SPDX-License-Identifier: MIT
//
// A minimal thread-safe logger. We deliberately avoid a heavy logging
// dependency; QMC runs primarily print progress / banners and do not
// produce massive log volumes.

#pragma once

#include <chrono>
#include <cstdio>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>

namespace qmc {

enum class LogLevel { Debug = 0, Info = 1, Warn = 2, Error = 3 };

class Log {
public:
    static Log& get() {
        static Log instance;
        return instance;
    }

    void set_level(LogLevel l) { level_ = l; }
    LogLevel level() const { return level_; }

    template <class... Args>
    void debug(Args&&... a) { write(LogLevel::Debug, std::forward<Args>(a)...); }
    template <class... Args>
    void info (Args&&... a) { write(LogLevel::Info,  std::forward<Args>(a)...); }
    template <class... Args>
    void warn (Args&&... a) { write(LogLevel::Warn,  std::forward<Args>(a)...); }
    template <class... Args>
    void error(Args&&... a) { write(LogLevel::Error, std::forward<Args>(a)...); }

private:
    LogLevel   level_ = LogLevel::Info;
    std::mutex mu_;

    static const char* tag(LogLevel l) {
        switch (l) {
            case LogLevel::Debug: return "DEBUG";
            case LogLevel::Info:  return "INFO ";
            case LogLevel::Warn:  return "WARN ";
            case LogLevel::Error: return "ERROR";
        }
        return "?    ";
    }

    template <class... Args>
    void write(LogLevel l, Args&&... a) {
        if (l < level_) return;
        std::ostringstream os;
        (os << ... << a);
        std::lock_guard<std::mutex> lock(mu_);
        std::fprintf(stderr, "[qmc:%s] %s\n", tag(l), os.str().c_str());
    }
};

} // namespace qmc
