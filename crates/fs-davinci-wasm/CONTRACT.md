# fs-davinci-wasm — CONTRACT

Layer: L6 browser boundary. This standalone workspace follows the narrow
patent-binding pattern while leaving generic joint semantics in `fs-mbd`.

## Purpose and layer

The boundary between the US 6,331,181 museum exhibit and
`fs_mbd::davinci::step_davinci_topology`. It exports the source-bounded
revolute/prismatic topology of Figs. 2 and 2A plus the claimed compatibility
identifier predicate. It does not invent dimensions, motor power, torque,
friction, accuracy, attenuation, clinical performance, or a universal
commercial robot specification.

## Public types and semantics

| Entry | Kind | Contract |
|---|---|---|
| `davinci_topology_step(...) -> String` | wasm + native | Admits finite radians, normalized insertion in `[-1,1]`, and compatibility-identifier presence; returns six generic joint axes/coordinates in `{"ok":{...}}`; every refusal uses `{"refusal":{...}}`. |

## Invariants

1. Five coordinates are generic revolute joints and insertion is one generic
   prismatic joint along the normalized tool axis.
2. Compatibility is reported as the patent's interface predicate. It does not
   silently rewrite joint state or claim a motor safety interlock.
3. Insertion is dimensionless display topology, not metres.
4. The boundary returns no source-unavailable force, speed, power, error,
   precision, or clinical metric.

## Error model

`input-outside-domain` refuses non-finite angles or insertion outside
`[-1,1]`. `multibody-refusal` preserves an unexpected generic joint-owner
refusal. Refusals carry a concise explanation and repair.

## Determinism class

Pure, synchronous, byte-identical output for identical inputs on the same
build. The path has no clock, randomness, host identity, or hidden state.

## Cancellation behavior

The bounded single-state query is synchronous and has no cancellation scope.

## Unsafe boundary

None. The crate forbids unsafe Rust.

## Feature flags

None. The wasm binding is selected only by target architecture.

## Conformance tests

Native tests cover joint count/axes, compatibility predicate semantics, and
out-of-domain refusal. The `fs-mbd` owner separately checks the generic joint
composition.

## No-claim boundaries

The connected exhibit geometry and host-side cup contact are pedagogical
models. This browser seam is not a calibrated commercial manipulator, safety
controller, or clinical device model.
