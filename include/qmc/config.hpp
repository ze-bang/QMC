// SPDX-License-Identifier: MIT
//
// A tiny dependency-free `key = value` configuration parser. The format
// is intentionally minimal so we don't need to vendor a JSON / TOML
// library:
//
//     # comments start with '#' and run to end of line
//     [section]                    # optional section headers (ignored)
//     key      = value             # whitespace around '=' is trimmed
//     list_key = 1, 2, 3           # comma separated lists
//
// Values are stored verbatim and converted on access.

#pragma once

#include <algorithm>
#include <cctype>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace qmc {

class Config {
public:
    Config() = default;

    static Config from_file(const std::string& path) {
        std::ifstream in(path);
        if (!in) {
            throw std::runtime_error("Config: cannot open file '" + path + "'");
        }
        std::stringstream ss;
        ss << in.rdbuf();
        Config c;
        c.parse(ss.str(), path);
        return c;
    }

    static Config from_string(const std::string& text) {
        Config c;
        c.parse(text, "<string>");
        return c;
    }

    bool has(const std::string& key) const { return kv_.count(key) != 0; }

    const std::string& raw(const std::string& key) const {
        auto it = kv_.find(key);
        if (it == kv_.end()) {
            throw std::out_of_range("Config: missing key '" + key + "'");
        }
        return it->second;
    }

    template <class T>
    T get(const std::string& key) const {
        return convert<T>(raw(key), key);
    }

    template <class T>
    T get_or(const std::string& key, T fallback) const {
        auto it = kv_.find(key);
        if (it == kv_.end()) return fallback;
        return convert<T>(it->second, key);
    }

    template <class T>
    std::vector<T> get_list(const std::string& key) const {
        const std::string& s = raw(key);
        std::vector<T> out;
        std::string token;
        std::stringstream ss(s);
        while (std::getline(ss, token, ',')) {
            trim(token);
            if (!token.empty()) out.push_back(convert<T>(token, key));
        }
        return out;
    }

    // Set programmatically (handy in tests).
    void set(std::string key, std::string value) { kv_[std::move(key)] = std::move(value); }

    // Iterate all key/value pairs (order is unspecified).
    auto begin() const { return kv_.begin(); }
    auto end()   const { return kv_.end();   }

private:
    std::unordered_map<std::string, std::string> kv_;

    static void trim(std::string& s) {
        auto issp = [](unsigned char c) { return std::isspace(c) != 0; };
        while (!s.empty() && issp(s.front())) s.erase(s.begin());
        while (!s.empty() && issp(s.back()))  s.pop_back();
    }

    void parse(const std::string& text, const std::string& source) {
        std::stringstream ss(text);
        std::string line;
        std::size_t lineno = 0;
        while (std::getline(ss, line)) {
            ++lineno;
            // Strip comments.
            const auto hash = line.find('#');
            if (hash != std::string::npos) line.erase(hash);
            trim(line);
            if (line.empty()) continue;
            if (line.front() == '[' && line.back() == ']') continue; // section header
            const auto eq = line.find('=');
            if (eq == std::string::npos) {
                throw std::runtime_error(source + ":" + std::to_string(lineno) +
                                         ": expected 'key = value'");
            }
            std::string key   = line.substr(0, eq);
            std::string value = line.substr(eq + 1);
            trim(key);
            trim(value);
            // Strip surrounding quotes if present.
            if (value.size() >= 2 &&
                ((value.front() == '"' && value.back() == '"') ||
                 (value.front() == '\'' && value.back() == '\''))) {
                value = value.substr(1, value.size() - 2);
            }
            kv_[std::move(key)] = std::move(value);
        }
    }

    template <class T>
    static T convert(const std::string& s, const std::string& key) {
        std::stringstream ss(s);
        if constexpr (std::is_same_v<T, bool>) {
            std::string t = s;
            std::transform(t.begin(), t.end(), t.begin(),
                           [](unsigned char c) { return std::tolower(c); });
            if (t == "true" || t == "yes" || t == "on" || t == "1") return true;
            if (t == "false"|| t == "no"  || t == "off"|| t == "0") return false;
            throw std::runtime_error("Config: cannot parse boolean for key '" + key + "': " + s);
        } else if constexpr (std::is_same_v<T, std::string>) {
            return s;
        } else {
            T value{};
            ss >> value;
            if (!ss) {
                throw std::runtime_error(
                    "Config: cannot parse value for key '" + key + "': " + s);
            }
            return value;
        }
    }
};

} // namespace qmc
