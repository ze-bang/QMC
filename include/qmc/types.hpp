// SPDX-License-Identifier: MIT
//
// Common scalar / index types used throughout the library.
//
// Keeping a single header for typedefs makes it trivial to switch the
// floating-point precision (e.g. to long double for very low-temperature
// runs) or to widen the integer index types for very large lattices.

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

namespace qmc {

using Real    = double;
using Site    = std::int32_t;     // index of a lattice site
using Bond    = std::int32_t;     // index of a bond
using OpCode  = std::int32_t;     // packed (bond, type) for one SSE operator
using Length  = std::int32_t;     // operator string position / length

inline constexpr OpCode kIdentity = -1;

// Pack / unpack helpers for an SSE operator. Type 0 = diagonal,
// type 1 = off-diagonal. Identity is represented by kIdentity (= -1).
inline constexpr OpCode pack_op(Bond b, int type) noexcept {
    return static_cast<OpCode>((b << 1) | (type & 1));
}
inline constexpr Bond op_bond(OpCode op) noexcept { return op >> 1; }
inline constexpr int  op_type(OpCode op) noexcept { return op & 1; }
inline constexpr bool op_is_identity(OpCode op) noexcept { return op == kIdentity; }
inline constexpr bool op_is_diagonal(OpCode op) noexcept {
    return op != kIdentity && (op & 1) == 0;
}
inline constexpr bool op_is_offdiag (OpCode op) noexcept {
    return op != kIdentity && (op & 1) == 1;
}

// Convenience: a contiguous spin configuration (±1).
using SpinConfig = std::vector<std::int8_t>;

} // namespace qmc
