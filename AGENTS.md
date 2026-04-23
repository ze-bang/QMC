# Agent guidelines for qmc_sse

Conventions to follow when modifying this codebase:

- **Header-only.** The library is intentionally header-only (everything
  in `include/qmc/*.hpp`). Do not introduce a `src/` layer unless a
  compile-time gain motivates it.
- **No external dependencies.** No FetchContent, no `find_package` for
  third-party libraries (OpenMP is the only exception and is optional).
- **C++20.** Use modern features liberally (concepts, `if constexpr`,
  ranges) but stay within what GCC 11 / Clang 14 / MSVC 19.30 support.
- **Style.** Anonymous namespaces or `static` for file-local helpers.
  Trailing-comma-friendly init lists, 4-space indent, lines <= 100 cols.
  Avoid `using namespace std;` in headers.
- **Performance.** Prefer flat `std::vector` storage to nested
  containers; avoid allocations inside Monte Carlo inner loops.
- **Reproducibility.** All randomness goes through `qmc::Pcg32` seeded
  from the user-controlled `seed` config key. Each replica must use a
  deterministic, distinct seed derived from `seed + r * stride`.
- **Tests.** Add a focused test for any new physics estimator. Physics
  tests should compare against an analytic limit / ED with a tolerance
  of a few standard errors.
- **Algorithm correctness.** When in doubt, cite the reference (paper
  + section / equation) in a comment near the code.
