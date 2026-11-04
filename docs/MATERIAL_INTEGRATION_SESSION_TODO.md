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
- [x] Verify all three code commits are on main and their changed paths are scoped.

## 2. Material heat storage in theta integration and its adjoint

Related plan work: MR19, `frankensim-material-reality-ui2if.3.1`;
MR04, `frankensim-material-reality-ui2if.1.4`.

- [x] Generalize the existing capacity assignment for both thermal history and theta marching.
- [x] Reuse the existing material-card density/Cp resolution and lumped capacity matrix.
- [x] Accept uniform, per-element, and material-ID capacity without an implicit fallback.
- [x] Carry material source domains through initial, prescribed, and accepted temperatures.
- [x] Use the same prepared capacity in forward and exact source-scale adjoint paths.
- [x] Assemble the sourced adjoint at an admitted physical temperature.
- [x] Differentiate the volume source without scaling fixed Robin/Neumann loads.
- [x] Add an independent material-dependent diffusion/mode and gradient regression.
- [x] Add heterogeneous heating and source-support refusal regressions.
- [x] Add fixed-boundary, zero-volume-source derivative regressions.
- [x] Review public exports, source admission, units, and cancellation boundaries.
- [x] Commit the implementation, focused tests, and bounded contract changes.
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
- [x] Add uniform mechanical-forcing regressions and review source-free branch equivalence.
- [x] Review source units, endpoint quadrature, field literals, and cumulative semantics.
- [x] Commit the implementation, focused tests, and bounded contract changes.
- [ ] Run the focused native NSF tests when execution is available.

## 5. Beads and native verification: explicit outstanding blockers

- [x] Locate the canonical `.beads/issues.jsonl` and its current Git blob.
- [x] Try supported GitHub file, line-range, raw-file, and blob reads.
- [x] Read the current Beads and create precise in-progress child tasks through `br`.
- [x] Record the exact committed scope and verification limitations through `br update --append-notes`.
- [x] Keep all three tasks in progress because native regression acceptance is still pending.
- [x] Preserve existing assignees and parent task states; isolate the three new records for publication.
- [x] Run `br --no-db sync --flush-only` after real task updates (success).
- [x] Run available rustfmt syntax checks and focused whitespace checks.
- [ ] Run the pinned-toolchain native Rust tests after dependency/disk blockers are resolved.

Execution was initially unavailable. The recovered environment provides an
existing nightly-2026-08-31 formatter, used as an explicitly different-toolchain
syntax/formatting check. The repository pins nightly-2026-07-06. The native
constellation preflight reports all seven required sibling checkouts missing;
shared-disk exhaustion prevents a responsible full bootstrap. No native Rust
test is reported as executed.

### Committed implementations

- Flow source-energy accounting: `215a11c0fc4576b3668f990375f3a2d9340957bc`.
- Native sourced Maxwell walls: `7161a318a6056e2ff99409415ff293d83bac6d7b`.
- Shared theta capacity and adjoint: `7fb1ca17948951e325b6a856d3f581a34f6fe02c`.

### Tracker state and publication blocker

The recovered checkout and actual `br 0.7.6 --no-db` created these scoped tasks
and appended their exact implementation commits, test scope and pending
verification. All remain `in_progress`; broad parents were not closed:

- Theta capacity/adjoint: `frankensim-material-reality-ui2if.3.1.2`.
- Native Maxwell wall: `frankensim-material-reality-ui2if.5.3.1`.
- NSF source energy: `frankensim-material-reality-ui2if.1.6.1`.

The native tracker contains the updates, and `br --no-db sync --flush-only`
completed successfully. Remote publication is still blocked: the canonical
export is approximately 41 MB, while the GitHub connector rejects request
bodies above 16 MiB. The native Git dry run also refused because no usable
password/credential was available. No remote tracker update is claimed.

`docs/MATERIAL_INTEGRATION_BEADS.patch` preserves the three exact `br`-authored
records, including their notes, as an unapplied publication patch. Its baseline
is Beads blob `03501e63c9c7a51581e3c27b66270510b11c1fa9`; applying only these
additions produces blob `ac58a8b7b8c667f80a687224b7858cebc499889e`. Existing
records are preserved byte for byte. A Git apply dry run against the exact
baseline index passed; the canonical working tracker was not modified by that check. Incidental normalization of old close
timestamps/comment ordering by the local `br` import is excluded from the patch.

- [x] Preserve the exact task updates in a small committed patch.
- [ ] Apply/publish the tracker additions from an authenticated native Git host,
  after checking the current canonical export for concurrent changes.
- [ ] Execute the focused native regression commands below with the pinned
  toolchain and constellation checkouts.
- [ ] Close the scoped Beads only after their regression acceptance passes.

### Focused native checks still pending

```sh
cargo test -p fs-conduction --test material_heat_capacity --test transient --test thermal_history_capacity
cargo test -p fs-duct --lib
cargo test -p fs-couple --test material_duct_wall
cargo test -p fs-lbm --test nsf3_thermochemical
```

The existing equations were reviewed, including independent scalar/matrix
checks of thermal mode/heating derivatives and acoustic passivity. These do
not substitute for compiling and executing the Rust regressions. No phase
transition, nonlinear caloric sensitivity, broadband material qualification,
or full-project material universality is claimed by these bounded changes.
