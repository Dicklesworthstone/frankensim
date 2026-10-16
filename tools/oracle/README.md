# FrankenSim high-precision oracle

This nested Cargo workspace is the external development oracle required by
Bead `frankensim-extreal-program-f85xj.3.1`. It audits two claims that cannot be
established by the production implementations themselves:

- every declared `fs-math` ULP budget for `exp`, `expm1`, `ln`, `sin`, `cos`,
  `tan`, `atan2`, `erf`, `pow`, `sqrt`, and `tanh` against correctly rounded
  MPFR results at 256-bit precision;
- `fs-ivl::Interval` arithmetic, negation, absolute value, hull,
  intersection, certified width, and elementary enclosures against exact-f64
  endpoint and interior witnesses evaluated by the same independent
  high-precision engine.

## Isolation decision

The oracle uses `rug` as a safe Rust interface to GMP/MPFR. That dependency is
intentionally foreign to FrankenSim's production dependency policy, so this
directory declares its own `[workspace]`, owns a separate reviewed
`Cargo.lock`, and is absent from the root workspace members. Neither `fs-math`
nor `fs-ivl` depends on the oracle. Root `xtask check-deps` must continue to see
the production graph exactly as before.

This is a comparison harness, not an admitted runtime authority. Its DSR row
binds the source snapshot, lockfile, precision, deterministic sample count,
per-function histogram, maximum observed ULP error, and argmax input. The
finite sample family does not prove every binary64 input, MPFR does not
authenticate the surrounding run, and a green row does not replace the
cross-ISA determinism or interval-law batteries.

## Running

```bash
cargo test --locked --manifest-path tools/oracle/Cargo.toml
cargo run --locked --manifest-path tools/oracle/Cargo.toml
```

The binary emits deterministic JSON Lines and exits nonzero if a declared ULP
budget is exceeded, a special-value policy disagrees with MPFR, or an interval
fails to contain an exact high-precision point evaluation. `samples` defaults
to 4096 per audited family and may be increased explicitly:

```bash
cargo run --locked --manifest-path tools/oracle/Cargo.toml -- \
  --samples 16384 --precision-bits 320
```

Precisions below 200 bits and zero sample counts are refused before any audit
row is emitted.
