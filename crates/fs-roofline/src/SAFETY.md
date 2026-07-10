# SAFETY: fs-roofline/axes.rs FMA probe capsule (bead xlvx)

## What the unsafe does
One `#[target_feature(enable = "avx,fma")]` function whose body is the
same pure-safe FMA chain as the portable twin — no pointers, no
intrinsics, only `f64::mul_add` over a caller-owned `[f64; LANES]`.

## Why it exists
Baseline x86-64 has no compile-time FMA: `mul_add` lowers to a libm
CALL and the measured "peak FLOPs" axis read 1.0 GFLOP/s on a
Threadripper (call overhead, not FMA throughput), inflating every
downstream attainment (GEMM read 40×). The target_feature body lets
LLVM emit real vfmadd for the probe without enabling FMA globally.

## Why it is sound
The dispatcher verifies `avx` + `fma` via
`is_x86_feature_detected!` immediately before the only call site. The
body performs no memory access beyond the `&mut` array the caller
lends. On non-x86 the capsule does not exist.

## Compensating checks
`fma_bench_reports_positive_throughput` guards the axis; the fs-la
perf lane's attainment gate (≤ ~1.0 by construction against an honest
axis) is the end-to-end trip-wire that caught the original bug.
