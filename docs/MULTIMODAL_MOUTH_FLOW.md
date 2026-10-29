# Mouth flow from the original multimodal bore

`fs_duct::modal::outlet` exposes `mouth_response`, `mouth_response_tabulated`,
`mouth_response_load` and `mouth_response_flanged`. Each returns the original
input `ModalResponse` and the complete row-major map from inlet modal volume
flow to mouth modal volume flow. Mode zero is total volume flow; higher modes
retain the same wall-normalized radial basis as the original modal owner.

The calculation runs inside the existing backward reflection recurrence.
At the terminal plane it retains `H = Zc^-1 (I-R)`. Propagating toward the inlet
updates `H <- H E` with the original decaying `E = diag(exp(i k L))`. Junction
pressure and flow continuity give the downstream forward amplitudes from the
upstream ones; `H` composes that map. It never evaluates a growing evanescent
exponential, inverts the projection matrix, or solves a substitute scalar bore.
Only one additional N-by-N map is retained, not a station-by-station field.
The impedance-only entry points disable this extra calculation and preserve
their original arithmetic and public result type.

`MouthResponse::per_incident_plane_pressure(reference)` gives every mouth-mode
flow per unit incident plane pressure at the specified positive real reference.
It imposes zero higher-mode INLET flow, matching the Z00 closure used by
`MmLineBank`. This is different from driving the full reflection matrix with
zero higher-mode incident waves. The positive-frequency table path validates
its original mouth plane, radius and gas; a separate explicit-load evaluation
can use declared static impedance without extending the table below its band.

The flanged entry point retains the complete original exterior radiation matrix
and its interior/exterior-gas admission. The other entry points retain their
existing matched higher-mode mouth approximation. Cone endpoints are still the
original midpoint-radius staircase image: returned endpoint radii identify the
actual discrete stations, and convergence still requires refining that image.
This change neither alters that geometry nor certifies arbitrary imported loads
as causal. Mouth flow is not a microphone signal or an exterior radiation model.

Five Rust regression groups cover analytic resistive cylinders, long evanescent
spans, scalar-chain degeneration in both junction directions, an independent
simultaneous all-wave system, flanged input/output power, and strict-table/DC
separation. Native Rust compilation and test execution remain unverified in the
implementation environment. An independent Python/SciPy derivation check passed
160 simultaneous-system comparisons (1-4 modes, both radius-step directions,
complex propagation and above/below-cutoff frequencies), with worst relative
flow error 2.03e-15, plus a long-evanescence case. That checks the equations, not
the compiled Rust implementation or experimental instrument accuracy.

Run the native regressions with:

```sh
cargo test -p fs-duct --lib modal::outlet::tests
```
