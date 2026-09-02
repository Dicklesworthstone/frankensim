# CONTRACT: fs-salisbury-wasm

## Purpose and layer

Layer: **L6 HELM/browser boundary**. This standalone crate exposes the
source-bounded US 4,921,293 hand composition from `fs-mbd` to a browser. The
generic owner remains `fs-mbd::salisbury`; this crate only admits scalar inputs
and serializes a stable JSON envelope.

## Public surface

`salisbury_hand_step(t1_n, t2_n, t3_n, t4_n, radius_scale_m,
first_idler_fixed)` returns a JSON string containing either:

- `ok`: the palm-rooted three-digit parent graph, nine-joint/12-cable topology,
  three generic axes, admitted SI tensions and radii, three source-law torques,
  Claim 1/2 predicates, and an explicit
  `historical_dynamics_available: false`; or
- `refusal`: a stable code, diagnostic, and repair.

## Invariants and determinism

- Cable tensions are finite and non-negative.
- The visitor-declared radius scale is finite and positive.
- The parent graph is exactly three palm-anchored serial three-joint chains.
- No partial or non-finite result is serialized.
- Successful output is a deterministic function of the admitted inputs on the
  same Rust/WASM toolchain.

## Cancellation and unsafe boundary

The computation is bounded and synchronous. It starts no tasks, performs no
I/O, and has no cancellation seam. Unsafe Rust is forbidden.

## No-claim boundary

This boundary does not supply historic dimensions, link geometry, masses,
inertias, damping, motor limits, cable material, friction, contact mechanics,
grasp force, force closure, stability, speed, or hardware validation. The
radius ratios are explicitly illustrative visitor-study inputs.

## Conformance

Unit tests cover successful topology/torque serialization, the Claim 2 boolean
probe, and typed refusal for invalid physical inputs. `fs-mbd` owns the G0
equation and joint-topology tests.
