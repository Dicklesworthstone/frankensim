# fs-daimler-wasm — CONTRACT

Layer: L6 browser boundary. This standalone workspace follows the narrow
patent-binding pattern while leaving the governing multibody semantics in the
generic `fs-mbd` owner.

## Purpose and layer

The boundary between the US 361,931 museum exhibit and
`fs_mbd::daimler::step_daimler_marine`. It exports the source's discrete
longitudinal shaft selections, ahead/astern contact topology, continuously
one-direction motor sign, and passive/optional-pump cooling paths. It does not
invent a shaft travel, friction coefficient, thrust, speed, cooling flow,
temperature, power, or efficiency.

## Public types and semantics

| Entry | Kind | Contract |
|---|---|---|
| `daimler_marine_step(shaft_selection, cooling_pump_enabled) -> String` | wasm + native | Admits exactly `-1` (astern), `0` (neutral), or `1` (ahead), composes the generic `fs-mbd` prismatic joint, and returns normalized shaft translation, its unit axis/DoF count, mutually exclusive drive contacts, motor/propeller direction signs, and cooling-path state inside `{"ok":{...}}`; every refusal uses `{"refusal":{...}}`. |

## Invariants

1. Ahead translates toward the motor (negative along the declared positive
   sternward axis), engages only `a/a²`, and permits the source-stated thrust
   contact predicate.
2. Astern translates away from the motor, opens `a/a²`, and engages only the
   `e¹/e²` with `a²/c` reverse path. Neutral opens both.
3. The motor direction sign is positive in all three states; only downstream
   contact topology changes propeller direction.
4. Enabling centrifugal pump `u` is additive and never erases the fore/aft
   outside-water pipe path printed in Claims 7–9.
5. The emitted shaft coordinate is normalized display topology, not metres.

## Error model

`input-outside-domain` refuses every shaft selector outside the three printed
reader states. `multibody-refusal` preserves an unexpected refusal from the
generic owner. Refusals carry a concise explanation and ranked repair.

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

Native tests cover ahead, neutral, astern, optional-pump composition, exact
mutual exclusion, and out-of-domain refusal. The `fs-mbd` owner tests separately
cover the generic joint and source-state contracts.

## No-claim boundaries

This renderer seam is not a calibrated clutch, propeller, cooling, or vessel
model. It returns no force, torque, coefficient, travel distance, rate, energy,
or performance quantity because US 361,931 supplies none of the required
inputs. Normalized visual distances remain presentation-owned.
