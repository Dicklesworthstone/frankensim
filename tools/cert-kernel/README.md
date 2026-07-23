# FrankenSim independent certified-arithmetic kernel

This nested Cargo workspace is the deliberately small N-version checker for
Bead `frankensim-extreal-program-f85xj.3.6`. It implements only interval
addition, subtraction, multiplication, division, square root, exponential, and
natural logarithm.

## Independence boundary

The kernel core has no dependency on `fs-ivl`, `fs-math`, MPFR, or a platform
elementary-function result:

- addition uses an exact `TwoSum` residual to choose the directed endpoint;
- multiplication, division, and square root compare exact binary64 dyadics in
  an `u128` significand representation;
- exponential uses power-of-two reduction, a positive Taylor series with a
  geometric tail bound, and interval squaring;
- logarithm uses exact binary range reduction and the positive
  `atanh` series with an explicit tail bound.

The feature-gated `crosscheck` module imports `fs-ivl` only after both kernels
have produced their answers. It checks nonempty intersections, exact
hand-derived references, and deterministic width ratios. A seeded mutation
shrinks the upper endpoint of `1 + 2^-53` by one ULP and proves the tripwire
detects the resulting false enclosure.

This is a diagnostic check, not a second production arithmetic authority.
Divergence is adjudicated against the exact-rational corpus and the isolated
MPFR lane; agreement cannot promote a certificate or evidence color.

## Running

```bash
cargo fmt --manifest-path tools/cert-kernel/Cargo.toml --check
cargo clippy --locked --manifest-path tools/cert-kernel/Cargo.toml \
  --all-targets -- -D warnings
cargo test --locked --manifest-path tools/cert-kernel/Cargo.toml
cargo run --quiet --locked --manifest-path tools/cert-kernel/Cargo.toml -- \
  --samples 4096
```

The binary emits deterministic JSON Lines. Each operation row includes exact
reference counts, compatibility failures, and width-ratio quantiles. It exits
nonzero on any exact-reference miss, non-overlap, or failed seeded tripwire.
