// SPDX-License-Identifier: MIT
#include <string>

#include "test_runner.hpp"

int main(int argc, char** argv) {
    std::string filter = (argc > 1) ? argv[1] : "";
    return qmc::test::run_all(filter);
}
