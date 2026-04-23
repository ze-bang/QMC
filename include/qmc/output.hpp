// SPDX-License-Identifier: MIT
//
// Minimal output helpers: timestamps, banner printing, and CSV writers
// for time series of measurements.

#pragma once

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "qmc/types.hpp"

namespace qmc {

inline std::string timestamp() {
    const auto now   = std::chrono::system_clock::now();
    const auto t     = std::chrono::system_clock::to_time_t(now);
    std::tm    tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

class CsvWriter {
public:
    CsvWriter(const std::string& path, std::vector<std::string> columns)
        : columns_(std::move(columns)) {
        out_.open(path);
        if (!out_) throw std::runtime_error("CsvWriter: cannot open " + path);
        for (std::size_t i = 0; i < columns_.size(); ++i) {
            if (i) out_ << ',';
            out_ << columns_[i];
        }
        out_ << '\n';
    }

    void write_row(const std::vector<Real>& row) {
        if (row.size() != columns_.size()) {
            throw std::runtime_error("CsvWriter: row width mismatch");
        }
        for (std::size_t i = 0; i < row.size(); ++i) {
            if (i) out_ << ',';
            out_ << std::scientific << std::setprecision(9) << row[i];
        }
        out_ << '\n';
    }

    void flush() { out_.flush(); }

private:
    std::vector<std::string> columns_;
    std::ofstream            out_;
};

} // namespace qmc
