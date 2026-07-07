# CONTRACT: fs-ladder

The fidelity-ladder registry (plan addendum, Proposal 3): the shared substrate
that the speculation proposer, the query planner, the Goodhart guard, and
tolerance allocation all walk.

## Purpose and layer

Layer L3 (FLUX-adjacent). Pure abstraction — NO numerical runtime dependencies
(std only). Each kernel declares an ordered ladder of fidelity rungs with
pluggable transfer operators between adjacent rungs.

## Public types and semantics

- `Rung { index, name, relative_cost, note }` — a fidelity level; index 0 is
  the coarsest/cheapest. Rungs are TOTALLY ORDERED per kernel.
- `Transfer` (trait) — `prolongate(coarse) -> fine` (rung k→k+1) and
  `restrict(fine) -> coarse` (rung k+1→k). Consumers supply the operator
  (fs-feec transfer, a correlation model, …).
- `Ladder` — builder `Ladder::new(kernel, base_name, base_cost, base_note)`
  then `.then(transfer, name, cost, note)` per finer rung (keeps
  `transfers.len() == rungs.len()-1` by construction). Query: `rung(k)`,
  `top()`, `bottom()`, `adjacent_rungs(k) -> AdjacentRungs { coarser, finer }`
  (empty in the correct direction at the ends), `prolongate(from, coarse)`,
  `restrict(from, fine)`.
- `LadderRegistry` — `register(ladder)`, `ladder(kernel)`, `kernels()`;
  `cht()` seeds the Proposal-7 conjugate-heat-transfer ladder
  (`correlation-Nu` → `RANS` → `LES`) so the ladder is immediately real.
- `Refine1d` — a concrete 1D coarsen/refine-by-2 demonstrator: prolongation is
  linear interpolation (`n → 2n-1`), restriction is injection (`2n-1 → n`), so
  `restrict ∘ prolongate = identity` and `prolongate ∘ restrict` is an
  idempotent projection.
- `LadderError` — `NoSuchRung` / `AtTop` (prolongate at the finest rung) /
  `AtBottom` (restrict at the coarsest rung) / `NoKernel`, each with a teaching
  `Display`.

## Invariants

- Rung indices are `0..len`, strictly increasing; `transfers[k]` bridges rung
  `k` and `k+1`.
- `adjacent_rungs` returns `coarser = None` at the bottom and `finer = None`
  at the top; out-of-range indices are structured errors, never panics.
- For a consistent transfer, `restrict(prolongate(x)) == x` (bit-exact for
  `Refine1d`).

## Error model

Structured `LadderError` values (a refusal that teaches), never panics, for
out-of-range rungs, prolongate-at-top, restrict-at-bottom, and unknown kernels.

## Determinism class

Fully deterministic: rung resolution and `Refine1d` transfers are pure
functions (no RNG, no I/O); outputs are bit-identical on replay — load-bearing
because transfer outputs feed verified-color certificates (G5).

## Cancellation behavior

None — all operations are bounded, synchronous pure functions (no `Cx`).

## Unsafe boundary

None. `#![deny(unsafe_code)]` via the workspace lint.

## Feature flags

None.

## Conformance tests

`tests/ladder.rs` (Proposal 3, 12 cases): total rung ordering; out-of-range
rung structured error; adjacency empty in the correct direction at the ends;
the G0 approximation property (`restrict∘prolongate = identity`,
`prolongate∘restrict` idempotent); prolongate-at-top / restrict-at-bottom
structured errors; registry register/resolve + unknown-kernel error; the CHT
ladder (`correlation-Nu`→`RANS`→`LES`); `Refine1d` degenerate lengths; G5
determinism.

## No-claim boundaries

- The registry owns rung DECLARATIONS + ordering + adjacency; the numerical
  `Transfer` between rungs is consumer-supplied. `Refine1d` is a 1D
  demonstrator, not a production multigrid transfer — real per-kernel operators
  (fs-feec prolongation/restriction, correlation↔RANS state maps) are injected.
- Rung `relative_cost` is an advisory hint; the query planner learns real costs
  from telemetry.
- The ladder does not run solves — it declares rungs and transfers states; the
  physics kernels execute the rungs.
