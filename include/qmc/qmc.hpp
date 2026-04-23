// SPDX-License-Identifier: MIT
//
// Umbrella header that pulls in the public API of the SSE QMC library.

#pragma once

#include "qmc/config.hpp"
#include "qmc/heisenberg.hpp"
#include "qmc/lattice.hpp"
#include "qmc/logging.hpp"
#include "qmc/measurements.hpp"
#include "qmc/observable.hpp"
#include "qmc/operator_string.hpp"
#include "qmc/output.hpp"
#include "qmc/rng.hpp"
#include "qmc/sse_engine.hpp"
#include "qmc/types.hpp"

namespace qmc {
inline constexpr const char* kVersion = "0.1.0";
}
