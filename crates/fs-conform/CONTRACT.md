# CONTRACT: fs-conform

The restriction-map plugin conformance SDK (plan addendum, Proposal 7):
certify chart-to-chart converters into tiers by the sheaf axioms. Owns risk R6.

## Purpose and layer

Layer L2 (restriction maps / Rep Router edges). No runtime dependencies — pure
Rust; `fs-propcheck` is test-only.

## Public types and semantics

- `Converter` (object-safe trait) — `id`, `source_dim`, `target_dim`,
  `apply` (source→target), `adjoint` (declared transpose), `declared_error`.
- `ManufacturedCase { input, exact_output }`; `Composition { after, direct,
  probes }` (a functoriality witness); `ConformanceSuite { adjoint_pairs,
  manufactured, composition, identity, tolerance }` (`identity` is an optional
  witness — probes on which a converter claiming to be the identity must act as
  one; `None` for non-identity converters).
- Axiom checks: `check_functoriality` (`after∘self == direct`), `check_identity`
  (`id(x) == x`), `check_adjoint` (`⟨Ax,y⟩ == ⟨x,Aᵀy⟩`),
  `check_tolerance_honesty` (measured ≤ declared) → `(honest, measured)`.
- `certify(&converter, &suite) -> ConformanceReport` — awards a `Tier`
  (`Rejected` / `Bronze` / `Silver` / `Gold`); reaches a tier above `Rejected`
  ONLY by passing every supplied axiom, with the level set by the (honestly met)
  declared error. `ConformanceReport::certified()`.

## Invariants

- A converter that fails ANY axiom (functoriality witness, adjoint, tolerance
  honesty) is `Rejected` — never certified.
- A converter that understates its error model FAILS tolerance honesty against
  manufactured solutions (a dishonest converter cannot buy a tier).
- Tier is monotone in declared-error tightness among conformant converters.
- The same `certify` path applies to first-party and third-party converters
  (R6: identical severity).

## Error model

Total functions returning booleans / a report (never panics); a rejection is a
report with `Tier::Rejected` and findings, not an error.

## Determinism class

Fully deterministic: all checks are pure functions of the converter + suite.

## Cancellation behavior

None (synchronous pure functions).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/conform.rs`: nine fixed Proposal 7 regression tests, two generated G0
laws, and one generated G3 conversion-path relation (512 cases each,
deterministic seeds and integrated shrinking). The fixed tests cover adjoint
consistency, tolerance honesty, correct and incorrect functoriality witnesses,
identity recognition and false-identity rejection, tiering, and uniform
first-party/third-party R6 severity. The generated G0 laws exercise exact
functor composition and identity action between the fixed pins. The G3 adopter
declares `restriction-map-direct-vs-composed` with
`conversion_path_independence` and applies it to real `Converter::apply`
routes (seed `0xC0F0_4A48_0003`, 512 cases, exact component tolerance): a
manufactured direct converter versus the staged `f(g(x))` route over nonzero,
exactly representable small-integer cases. The joint shrink surface retains the
matrix/probe input and route transform.

## No-claim boundaries

- Converters are modeled as finite-dimensional linear operators for the
  conformance harness; the SDK surface (the `Converter` trait) is what a real
  chart-to-chart trace/conversion operator implements.
- The generated G0 laws use small, exactly representable 2x2 integer matrices
  and probes. They do not certify arbitrary nonlinear converters or general
  floating-point associativity.
- The G3 path-independence adopter covers the same finite-dimensional linear
  fixture family. Exact component equality is justified by its bounded integer
  arithmetic; it does not extend that tolerance to nonlinear, approximate, or
  floating-roundoff-sensitive chart conversion paths.
- The suite is SUPPLIED here (probes, manufactured cases, composition witness);
  AUTO-GENERATING it from a chart pair's sheaf axioms is the generator's job (a
  downstream producer feeding this harness).
- Stamping the certified tier onto every ledger entry the converter touches is
  fs-ledger's integration; this crate produces the tier.
- Adjoint consistency is checked against the converter's OWN declared adjoint;
  cross-checking it against Proposal 1's ledger adjoint is the ledger's wiring.
