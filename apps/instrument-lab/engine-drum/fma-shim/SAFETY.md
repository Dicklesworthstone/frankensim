# fs-live-fma: unsafe capsule

## Why this exists

FrankenSim's deterministic elementary functions (`fs_math::det::{sin, cos, exp}`)
are Horner chains of `f64::mul_add`. The prepared exact-ZOH drum engine calls
them for every mode on every mechanical step. WebAssembly has no scalar fused
multiply-add instruction, so rustc lowers `mul_add` to a call to the software
`fma` of compiler-builtins. A Node CPU profile of the baseline build put
**79.8 % of the engine's time inside that one function**.

WebAssembly *relaxed-simd* has `f64x2.relaxed_madd`, which an engine may
implement fused or unfused. On hardware with FMA (every ARM64 CPU; x86-64 since
Haswell/Piledriver) V8, SpiderMonkey and JavaScriptCore emit the CPU's fused
instruction.

## What it does

`fma(a, b, c)` (the C symbol the compiler calls) is defined here:

1. On first use it probes the host once with an input whose fused and unfused
   results differ in every bit that matters:
   `(1 + 2^-30) * (1 + 2^-30) - (1 + 2^-29)` is exactly `2^-60` when fused and
   exactly `0` when the product is rounded first.
2. If the probe proves fusion, every later call is one `f64x2.relaxed_madd`
   on lane 0. A hardware FMA is correctly rounded by IEEE 754, so the result is
   bit-identical to the software routine it replaces and to native builds.
3. If the probe does not prove fusion, every call goes to `libm::fma`, the same
   correctly rounded software algorithm as the baseline build.

No numerical result of FrankenSim changes in either case. The engine reports
which path is active (`fma` in `info()`), and the worker falls back to the
baseline `.wasm` (built without this crate) when the browser cannot compile
relaxed-simd at all.

## The two unsafe items

- `#[unsafe(no_mangle)] extern "C" fn fma`: exports the libcall symbol. It
  must have exactly the C signature `double fma(double, double, double)`; it
  does. compiler-builtins defines the same symbol weakly, so this strong
  definition wins at link time.
- One call to a `#[target_feature(enable = "relaxed-simd")]` function. The
  module is only ever instantiated by engines that validated its relaxed-simd
  opcodes (validation is all-or-nothing for a WebAssembly module), so the
  feature is present whenever the code can run.

## Evidence

`engine-drum/tests` and the browser test page compare the hardware-FMA build
against the baseline build bit-for-bit on rendered audio (see the report).
