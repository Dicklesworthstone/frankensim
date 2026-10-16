# fs-tilelang-macros CONTRACT

## Purpose and layer

Layer: **UTIL** (proc-macro crate; zero dependencies — std
`proc_macro` only, the Franken-only law covers macro deps). The
`kernel!` function-like macro: parses the restricted kernel grammar,
performs static analysis, and generates the scalar/lane variants,
metadata, and twin tests that fs-tilelang's runtime types describe.

## Public types and semantics

`kernel! { name, reads, index_reads?, uparams?, params?, writes,
halo?, reduction, body }`:

- `reads` are `&[f64]`, `index_reads` are `&[u32]` (rewritten to
  usize via checked conversion), `uparams` are `usize` scalars,
  `params` are `f64` scalars, `writes` are `&mut [f64]`.
- Body grammar: `let` locals and assignments over bare buffer names
  (rewritten to indexed accesses at the implicit element index);
  `shift_add(buf, EXPR)` / `shift_sub(buf, EXPR)` stencil accesses;
  `gather(buf, EXPR)` computed-index reads; `acc = EXPR;` per-element
  reduction contribution (requires a declared reduction).
- Generated per kernel: `META` (flops counted from body operators,
  `mul_add` = 2; bytes = 8/f64 + 4/u32 buffer; literal halos recorded,
  dynamic halos as 0), `run_scalar`, `run` (single tier dispatch),
  `run_lanes::<LANES>`, and — when the macro can drive the kernel
  safely (no gather, no uparams) — a `#[cfg(test)]` twin-test module
  asserting bitwise tier equivalence at LANES ∈ {2, 4, 8} and repeat
  determinism.
- Length policy: writes anchor the common element count; bare-named
  and shift-target reads plus index reads must match (asserted at
  entry); gather-only buffers keep their own length.

## Invariants

- Expansion is a pure function of the input tokens (no environment,
  no randomness) — deterministic.
- Scalar and lane variants are generated from the SAME rewritten body
  string: per-element arithmetic cannot diverge between variants.
- Generated code compiles under the workspace clippy all+pedantic
  wall and `missing_docs`.
- Token rendering respects `Spacing::Joint` (`::`, `->` never split).

## Error model

Structured `compile_error!` with actionable messages — no silent
fallbacks: read/write aliasing and duplicate identifiers, reserved
names (`acc`, `gather`, `shift_add`, `shift_sub`), unsafe blocks,
user loops (`for`/`while`/`loop`), allocation attempts
(`vec!`/`Vec`/`Box`/`String`/`collect`/`push`), early `return`,
unknown declarations or reductions, `acc` without a reduction and
reductions without `acc`, writes never assigned, `shift`/`gather` on
undeclared buffers, malformed lists.

## Determinism class

Compile-time only; deterministic expansion per the invariant above.

## Cancellation behavior

Not applicable (runs inside rustc, single bounded pass).

## Unsafe boundary

None. `unsafe_code = "deny"`.

## Feature flags

None.

## Conformance tests

Exercised through fs-tilelang (a proc-macro crate cannot unit-test
its own expansion): `crates/fs-tilelang/tests/tilelang_battery.rs`
covers map/stencil/gather/reduction kernels end to end including the
auto-generated twin tests; `crates/fs-tilelang/tests/compile_fail.rs`
pins all ten diagnostic paths via the in-house offline harness.

## No-claim boundaries

- No span-precise diagnostics (string-based generation; errors point
  at the macro call).
- The flop counter is a token census (counts operator tokens and
  `mul_add`), not a dataflow analysis: dead code inflates it,
  function calls other than `mul_add` are uncounted. Good enough for
  roofline BINNING; not a cycle model. Bytes-per-element counts each
  declared buffer once (gather traffic is data-dependent and
  approximated the same way).
- One implicit element index; no 2D/3D index spaces, no cross-element
  communication, no writes at shifted positions (scatter), no
  non-f64 write buffers.
- Alias analysis is name-identity (buffers are distinct parameters;
  the CALLER's borrow checker enforces actual non-overlap of the
  slices passed in — two `&mut` can't alias in safe Rust).
