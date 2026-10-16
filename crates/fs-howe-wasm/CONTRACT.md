# fs-howe-wasm — CONTRACT

Layer: L6 browser boundary. This standalone workspace follows the narrow
patent-binding pattern while leaving generic joint semantics in `fs-mbd`.

## Purpose and layer

This crate binds the US 4,750 museum exhibit to
`fs_mbd::howe::step_howe_topology`. It exports the one-shaft, source-order
needle/loop/shuttle/feed topology and Claim 1 interlock predicate. It does not
invent overall dimensions, speed, mass, inertia, force, torque, friction,
thread strength, seam strength, power, or productivity.

## Public types and semantics

| Entry | Kind | Contract |
|---|---|---|
| `howe_topology_step(...) -> String` | wasm + native | Admits a finite main-shaft angle, normalized loop slack in `[0,1]`, and Claim 1 presence; returns generic joint axes, normalized coordinates, phase predicates, and the two printed local dimensions in `{"ok":{...}}`; every refusal uses `{"refusal":{...}}`. |

## Invariants

1. Seven scalar joint coordinates are constrained by the printed cams and
   linkages to one prescribed main-shaft drive coordinate.
2. The curved needle remains part of the vibrating arm, shuttle K remains on a
   prismatic axis in trough I, and baster plate H owns the cloth feed.
3. The 1/8-inch needle-eye offset and 3/4-inch baster-point pitch are the only
   source dimensions returned.
4. Normalized loop slack and a counterfactual Claim 1 offset are presentation
   controls, not historical measurements.

## Error model

`input-outside-domain` refuses non-finite input or loop slack outside `[0,1]`.
`multibody-refusal` preserves an unexpected generic joint-owner refusal. Each
refusal includes a concise repair.

## Determinism and cancellation

The bounded query is pure, synchronous, and byte-identical for identical input
on the same build. It has no clock, randomness, host identity, hidden state, or
cancellation scope.

## Unsafe boundary and features

None. The crate forbids unsafe Rust. The wasm binding is selected only by target
architecture and the crate has no optional features.

## Conformance and no-claim boundaries

Native tests cover axes, one-drive composition, source-order loop passage,
Claim 1 removal, printed dimensions, and invalid-input refusal. This is a
source-bounded kinematic exhibit, not a calibrated surviving-machine model or a
seam-strength, production-rate, or energy simulation.
