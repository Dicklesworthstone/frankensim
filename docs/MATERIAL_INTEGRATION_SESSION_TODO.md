# Material integration execution checklist

Starting revision: `dbb69a21fef6bfedb6febb47887c51efe5ed2bae`.

This checklist records the bounded implementation work selected for this session.
A checked implementation item means the source change was reviewed and committed;
it does not mean an unexecuted Rust test passed. The canonical task tracker remains
Beads. After execution access recovered, its live export was read through `br`.

## 1. Select useful work and preserve current main

- [x] Read the operating contract and relevant comprehensive-plan sections.
- [x] Read the Material Reality owners, physical requirements, and task mapping.
- [x] Compare main with the previous session's published revision.
- [x] Preserve the newer species-recirculation and gas source-domain changes.
- [x] Inspect and preserve concurrent magnetic, aeroacoustic, and GPU-renderer commits.
- [x] Inspect existing source refresh, field thermodynamics, and state identities.
- [x] Choose actual missing solver connections rather than duplicate those owners.
- [x] Keep independent thermal, acoustic, and flow edits in disjoint file sets.
- [x] Begin publishing coherent implementations directly to main with expected-head leases.
- [ ] Verify the final main revision and the complete set of changed paths.

## 2. Material heat storage in theta integration and its adjoint

Related plan work: MR19, `frankensim-material-reality-ui2if.3.1`;
MR04, `frankensim-material-reality-ui2if.1.4`.

- [ ] Generalize the existing capacity assignment for both thermal history and theta marching.
- [ ] Reuse the existing material-card density/Cp resolution and lumped capacity matrix.
- [ ] Accept uniform, per-element, and material-ID capacity without an implicit fallback.
- [ ] Carry material source domains through initial, prescribed, and accepted temperatures.
- [ ] Use the same prepared capacity in forward and exact source-scale adjoint paths.
- [ ] Assemble the sourced adjoint at an admitted physical temperature.
- [ ] Differentiate the volume source without scaling fixed Robin/Neumann loads.
- [ ] Add an independent material-dependent diffusion/mode and gradient regression.
- [ ] Add heterogeneous heating and source-support refusal regressions.
- [ ] Add fixed-boundary, zero-volume-source derivative regressions.
- [ ] Review public exports, source admission, units, and cancellation boundaries.
- [ ] Commit the implementation, focused tests, and bounded contract changes.
- [ ] Run the focused native conduction tests when execution is available.

## 3. Source-derived Maxwell relaxation in native duct dynamics

Related plan work: MR33, `frankensim-material-reality-ui2if.5.3`.

- [x] Add existing Maxwell branches to the native DG wall operator.
- [x] Preserve reciprocal pressure/wall coupling and every existing gas-memory block.
- [x] Include branch storage in wall observations and total native energy.
- [x] Include branch relaxation loss exactly once in dissipated power.
- [x] Keep ports, operator identity, and the existing DOF cap consistent.
- [x] Wire `LocalMaterialWall::bind_dg` to its retained Maxwell material spectrum.
- [x] Retain the source frequency band and current cylinder/gas admission limits.
- [x] Add an independent sealed-gas relaxing-mode regression.
- [x] Add zero-branch equivalence and combined-loss/port/state regressions.
- [x] Add a source-only relaxation-time change through the material adapter.
- [x] Review state indexing, branch units, source admission, and energy signs.
- [x] Commit the implementation, focused tests, and bounded contract changes.
- [ ] Run the focused native duct and coupling tests when execution is available.

## 4. External heating in the native compressible-flow energy account

Related plan work: MR06, `frankensim-material-reality-ui2if.1.6`;
MR07, `frankensim-material-reality-ui2if.1.7`.

- [x] Integrate exactly the source-energy quadrature used by the existing NSF step.
- [x] Expose signed supplied energy in the step and cumulative accounts.
- [x] Subtract supplied energy from the independent total-energy residual.
- [x] Exclude pure heat supply from the residual labelled viscous heating.
- [x] Distinguish unresolved internal exchange for manufactured mass/momentum sources.
- [x] Keep unresolved exchange out of the total-energy residual.
- [x] Preserve arithmetic and observations for the source-free path.
- [x] Retain atomic state/account rollback on material-support refusal.
- [x] Add constant heating, cooling, and time-varying source regressions.
- [x] Add uniform mechanical-forcing and source-free equivalence regressions.
- [x] Review source units, endpoint quadrature, field literals, and cumulative semantics.
- [x] Commit the implementation, focused tests, and bounded contract changes.
- [ ] Run the focused native NSF tests when execution is available.

## 5. Beads and native verification: explicit outstanding blockers

- [x] Locate the canonical `.beads/issues.jsonl` and its current Git blob.
- [x] Try supported GitHub file, line-range, raw-file, and blob reads.
- [x] Read the current Beads and create precise in-progress child tasks through `br`.
- [ ] Record the exact committed scope and verification limitations on those Beads.
- [ ] Close only tasks whose complete acceptance has actually been met.
- [ ] Preserve existing assignees and closed antecedents.
- [ ] Run `br sync --flush-only` after real task updates.
- [ ] Run applicable rustfmt and focused native Rust verification.

Execution was initially unavailable. The recovered environment provides an
existing nightly-2026-08-31 formatter, used as an explicitly different-toolchain
syntax/formatting check. The repository pins nightly-2026-07-06. The native
constellation preflight reports all seven required sibling checkouts missing;
shared-disk exhaustion prevents a responsible full bootstrap. No native Rust
test is reported as executed.

The 40,991,716-byte Beads export could not be fetched through the GitHub file
connector, but the recovered checkout and actual `br 0.7.6 --no-db` now work.
The existing broad parents retain their prerequisite state. Three scoped child
tasks were created in progress:

- Theta capacity/adjoint: `frankensim-material-reality-ui2if.3.1.2`.
- Native Maxwell wall: `frankensim-material-reality-ui2if.5.3.1`.
- NSF source energy: `frankensim-material-reality-ui2if.1.6.1`.

Their exact commit notes, final native-verification state, and exported tracker
publication are still being completed. No broader task is claimed closed.
