# CONTRACT: fs-crump-wasm

## Purpose and layer

Layer: **L6 HELM/browser boundary**. This standalone crate exposes a narrow,
source-bounded physical screen for the flexible-strand embodiment of US
5,121,329. It composes laws owned by generic FrankenSim crates:

- `fs-flux::capillary` owns the Newtonian circular-capillary pressure, wall
  shear-rate, and hydraulic-power relation.
- `fs-conduction::reduced_slab` owns the fixed-boundary first-mode cooling
  time and threshold crossing.

This crate owns only scalar admission and stable JSON serialization.

## Public surface

`crump_fdm_step(...)` returns a JSON string containing either:

- `ok`: SI capillary and first-mode thermal-screen outputs with the two owner
  names and their applicability boundaries; or
- `refusal`: a stable code, diagnostic, and repair.

## Invariants and determinism

- Geometry and material quantities are finite and strictly positive.
- Volumetric flow is finite and non-negative.
- Temperatures describe monotone cooling through the declared threshold.
- No partial or non-finite result is serialized.
- Successful output is a deterministic function of admitted inputs on the
  same Rust/WASM toolchain.

## Cancellation and unsafe boundary

The computation is bounded and synchronous. It starts no tasks, performs no
I/O, and has no cancellation seam. Unsafe Rust is forbidden.

## No-claim boundary

The capillary result assumes incompressible, Newtonian, fully developed,
laminar, no-slip flow in a straight rigid circular land. It does not model
polymer shear thinning, viscoelasticity, wall slip, entrance/contraction loss,
or free-surface deposition. The thermal result is a fixed-boundary,
one-dimensional, first-mode exponential screen. It does not infer convection,
contact resistance, phase change, crystallization, bonding, strength, or
historic machine performance. A glass-transition crossing is not labeled
solidification.

## Conformance

Unit tests cover generic-owner serialization, zero-flow admission, exact
thermal threshold crossing, and typed refusal for invalid capillary and
temperature inputs. The two generic crates own their respective G0 equations.
